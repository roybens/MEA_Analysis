#!/usr/bin/env python3
"""
emulator.py  —  Shared Memory Streamer Emulator
================================================
Reads recording from HDF5 file and emulates real-time streaming to shared memory.

This simulates the behavior of mxw_streamer.cpp by:
  1. Opening HDF5 file saved by recorder.py
  2. Creating POSIX shared memory segment /mxw_scope
  3. Writing header with metadata
  4. Streaming data to ring buffer at real-time speed (20 kHz)

Usage
-----
python emulator.py /tmp/mea_recording.h5 --speed 1.0

Arguments
---------
file              Input HDF5 file path (from recorder.py output)
--speed           Playback speed multiplier (1.0 = real-time, >1.0 = faster)
--channels        Channels to stream (auto-detected if omitted)
--loop            Loop playback when done (default: false)
"""

import argparse
import mmap
import os
import struct
import sys
import time
from pathlib import Path

import numpy as np

# Try h5py first, fallback to other methods
try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False
    print("[WARNING] h5py not available — trying binary fallback", file=sys.stderr)

# POSIX IPC for shared memory
try:
    import posix_ipc
    HAS_POSIX_IPC = True
except ImportError:
    HAS_POSIX_IPC = False
    print("[ERROR] posix_ipc not available. Install: pip install posix_ipc", file=sys.stderr)
    sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# Shared Memory Layout (must match oscilloscope.py)
# ─────────────────────────────────────────────────────────────────────────────

SHM_NAME         = "/mxw_scope"
SHM_MAGIC        = 0xC0FFEE42
MAX_CHANNELS     = 8
RING_FRAMES      = 40_000           # 2 seconds @ 20 kHz
SAMPLING_RATE    = 20_000           # Hz
DATA_OFFSET      = 4096             # Header is first 4096 bytes
DATA_SIZE        = RING_FRAMES * MAX_CHANNELS * 4  # float32
SHM_TOTAL        = DATA_OFFSET + DATA_SIZE

# Header struct format (must match oscilloscope.py HEADER_FMT)
HEADER_FMT       = "=II8HQQIF"
HEADER_SIZE      = struct.calcsize(HEADER_FMT)


def parse_args():
    p = argparse.ArgumentParser(description="Emulate mxw_streamer via shared memory")
    p.add_argument("file", type=str,
                   help="Input HDF5 file from recorder.py")
    p.add_argument("--speed", type=float, default=1.0,
                   help="Playback speed (1.0=real-time, default: 1.0)")
    p.add_argument("--channels", type=int, nargs="+", default=None,
                   help="Channels to stream (default: auto-detect)")
    p.add_argument("--loop", action="store_true",
                   help="Loop playback when done")
    return p.parse_args()


class SharedMemoryEmulator:
    """Manages POSIX shared memory ring buffer."""

    def __init__(self, n_channels, sampling_rate, channel_ids, lsb_uv=1.0):
        """
        Initialize shared memory segment.
        
        Args:
            n_channels: number of channels to stream
            sampling_rate: sample rate in Hz
            channel_ids: list of channel indices
            lsb_uv: LSB calibration in µV (default 1.0)
        """
        self.n_channels = n_channels
        self.sampling_rate = sampling_rate
        self.channel_ids = list(channel_ids)
        self.lsb_uv = lsb_uv

        # Clean up any existing shared memory
        try:
            posix_ipc.SharedMemory(SHM_NAME).unlink()
        except (ValueError, FileNotFoundError):
            pass

        # Create new shared memory segment
        print(f"[emulator] Creating shared memory {SHM_NAME} ({SHM_TOTAL} bytes)")
        self.shm = posix_ipc.SharedMemory(SHM_NAME, posix_ipc.O_CREAT | posix_ipc.O_EXCL,
                                          size=SHM_TOTAL)
        self.buf = mmap.mmap(self.shm.fd, SHM_TOTAL)
        os.close(self.shm.fd)  # Close fd, mmap holds the reference

        # Initialize ring buffer to zeros
        self.buf.seek(0)
        self.buf.write(b'\x00' * SHM_TOTAL)
        self.buf.flush()

        # Write header
        self._write_header()
        self.write_pos = 0
        self.frame_number = 0

    def _write_header(self):
        """Pack and write header to shared memory."""
        # Ensure channel_ids list is exactly 8 elements (padded with zeros)
        ch_ids = self.channel_ids + [0] * (8 - len(self.channel_ids))
        ch_ids = ch_ids[:8]

        header_tuple = (
            SHM_MAGIC,                  # magic
            self.n_channels,            # n_channels
            *ch_ids,                    # channel_ids[8]
            0,                          # write_pos (will be updated)
            0,                          # frame_number (will be updated)
            self.sampling_rate,         # sample_rate
            self.lsb_uv,                # lsb_uv
        )
        header_bytes = struct.pack(HEADER_FMT, *header_tuple)
        self.buf.seek(0)
        self.buf.write(header_bytes)
        # Pad rest of header to DATA_OFFSET
        self.buf.write(b'\x00' * (DATA_OFFSET - len(header_bytes)))
        self.buf.flush()

    def write_frame(self, data):
        """
        Write one frame (all channels) to ring buffer.
        
        Args:
            data: 1D numpy array of length n_channels (float32)
        """
        # Determine write position in ring buffer
        write_slot = self.write_pos % RING_FRAMES
        offset_bytes = DATA_OFFSET + (write_slot * MAX_CHANNELS * 4)

        # Create frame buffer: pad data to MAX_CHANNELS
        frame = np.zeros(MAX_CHANNELS, dtype=np.float32)
        frame[:self.n_channels] = data[:self.n_channels]

        # Write to ring buffer
        self.buf.seek(offset_bytes)
        self.buf.write(frame.tobytes())
        self.buf.flush()

        # Update write_pos and frame_number in header
        self.write_pos += 1
        self.frame_number += 1
        pos_bytes = struct.pack("=Q", self.write_pos)
        frame_bytes = struct.pack("=Q", self.frame_number)
        
        self.buf.seek(struct.calcsize("=II8H"))  # Offset to write_pos field
        self.buf.write(pos_bytes)
        self.buf.write(frame_bytes)
        self.buf.flush()

    def cleanup(self):
        """Close and unlink shared memory."""
        try:
            self.buf.close()
            self.shm.unlink()
            print(f"[emulator] Cleaned up shared memory")
        except Exception as e:
            print(f"[WARNING] Cleanup error: {e}", file=sys.stderr)


class HDF5Reader:
    """Read data from HDF5 file."""

    def __init__(self, filepath):
        if not HAS_H5PY:
            raise ImportError("h5py required for HDF5 files")
        
        self.file = h5py.File(filepath, 'r')
        self._discover_structure()

    def _discover_structure(self):
        """Auto-detect HDF5 structure and dataset names."""
        print(f"[emulator] HDF5 structure:")
        
        def print_structure(name, obj):
            print(f"  {name}: {type(obj).__name__}", end="")
            if isinstance(obj, h5py.Dataset):
                print(f" shape={obj.shape} dtype={obj.dtype}", end="")
            print()
        
        self.file.visititems(print_structure)

        # Try common dataset names
        self.data_key = None
        for candidate in ["data", "recording", "raw", "signals", "electrode_data"]:
            if candidate in self.file:
                self.data_key = candidate
                break
        
        if not self.data_key:
            # Use first dataset found
            keys = [k for k in self.file.keys() if isinstance(self.file[k], h5py.Dataset)]
            if keys:
                self.data_key = keys[0]
        
        if not self.data_key:
            raise ValueError("No datasets found in HDF5 file")
        
        dset = self.file[self.data_key]
        print(f"[emulator] Using dataset '{self.data_key}' shape={dset.shape} dtype={dset.dtype}")

        # Determine data shape: assume (time, channels) or (channels, time)
        if dset.shape[0] > dset.shape[1]:
            # First dim larger → likely time
            self.time_axis = 0
            self.n_samples, self.n_channels = dset.shape
        else:
            # Second dim larger → likely time
            self.time_axis = 1
            self.n_channels, self.n_samples = dset.shape
        
        print(f"[emulator] Data shape: {self.n_samples} samples × {self.n_channels} channels")

    def read_frame(self, sample_index):
        """
        Read one frame (all channels) at given sample index.
        
        Returns:
            1D numpy array of length n_channels (float32)
        """
        dset = self.file[self.data_key]
        if self.time_axis == 0:
            frame = dset[sample_index, :].astype(np.float32)
        else:
            frame = dset[:, sample_index].astype(np.float32)
        return frame

    def get_duration(self, sample_rate=20_000):
        """Get recording duration in seconds."""
        return self.n_samples / sample_rate

    def close(self):
        self.file.close()


def main():
    args = parse_args()

    if not os.path.exists(args.file):
        print(f"[ERROR] File not found: {args.file}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("  Shared Memory Emulator (mxw_streamer replacement)")
    print("=" * 60)
    print(f"  Input file   : {args.file}")
    print(f"  Playback speed: {args.speed}×")
    print(f"  Loop         : {args.loop}")
    print("=" * 60)
    print()

    try:
        # Open HDF5 file
        print(f"[emulator] Opening {args.file}...")
        reader = HDF5Reader(args.file)

        n_channels = args.channels[0] if args.channels else reader.n_channels
        if args.channels:
            channel_ids = args.channels[:8]
        else:
            channel_ids = list(range(min(reader.n_channels, MAX_CHANNELS)))

        # Initialize shared memory
        shm = SharedMemoryEmulator(
            n_channels=n_channels,
            sampling_rate=SAMPLING_RATE,
            channel_ids=channel_ids,
            lsb_uv=1.0
        )

        print()
        print(f"[emulator] Starting playback...")
        print(f"  Channels    : {channel_ids}")
        print(f"  Duration    : {reader.get_duration():.1f} seconds")
        print(f"  Sample rate : {SAMPLING_RATE} Hz")
        print()

        # Calculate sleep time between frames (emulate 20 kHz)
        frame_duration = 1.0 / SAMPLING_RATE  # 50 µs per frame
        sleep_time = frame_duration / args.speed

        try:
            loop_count = 0
            while True:
                loop_count += 1
                print(f"[emulator] Playback loop {loop_count}")

                start_time = time.time()
                for sample_idx in range(reader.n_samples):
                    # Read frame from file
                    frame_data = reader.read_frame(sample_idx)

                    # Write to shared memory
                    shm.write_frame(frame_data)

                    # Emulate real-time timing
                    elapsed = time.time() - start_time
                    ideal_elapsed = (sample_idx + 1) * sleep_time
                    sleep_ms = (ideal_elapsed - elapsed) * 1000
                    if sleep_ms > 0:
                        time.sleep(sleep_ms / 1000.0)

                    # Progress indicator
                    if (sample_idx + 1) % (SAMPLING_RATE // 2) == 0:  # Every 0.5 sec
                        elapsed = time.time() - start_time
                        print(f"  {sample_idx + 1}/{reader.n_samples} "
                              f"({elapsed:.1f}s elapsed)", end="\r")

                print()
                print(f"[emulator] ✓ Loop {loop_count} complete")

                if not args.loop:
                    break

                print(f"[emulator] Looping... (Ctrl-C to exit)")
                time.sleep(2)

        except KeyboardInterrupt:
            print("\n[emulator] Interrupted by user")

        finally:
            reader.close()
            shm.cleanup()

        print()
        print(f"[emulator] Emulation finished. Shared memory still available for oscilloscope.py")
        print(f"           Press Ctrl-C in oscilloscope to exit.")

    except Exception as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
