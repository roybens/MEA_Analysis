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

# Try h5py for .h5 files
try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

# Try numpy for reading .npy or raw binary files  
try:
    import scipy.io as sio
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

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
HEADER_FMT       = "=IIHHHHHHHHQQIf"
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
            shm_existing = posix_ipc.SharedMemory(SHM_NAME)
            try:
                shm_existing.unlink()
                print(f"[emulator] Cleaned up existing shared memory segment")
            except (ValueError, OSError, posix_ipc.PermissionsError) as e:
                print(f"[WARNING] Could not unlink existing segment: {e}")
                print(f"[WARNING] Will try to reuse it...")
        except (ValueError, FileNotFoundError, posix_ipc.PermissionsError, posix_ipc.ExistentialError):
            # No existing segment or can't access it, continue
            pass

        # Create new shared memory segment (or attach to existing one)
        print(f"[emulator] Creating/attaching shared memory {SHM_NAME} ({SHM_TOTAL} bytes)")
        try:
            # Try to create new segment
            self.shm = posix_ipc.SharedMemory(SHM_NAME, posix_ipc.O_CREAT | posix_ipc.O_EXCL,
                                              size=SHM_TOTAL)
            print(f"[emulator] Created new shared memory segment")
            created_new = True
        except (FileExistsError, posix_ipc.ExistentialError):
            # Segment exists, try to attach to it
            print(f"[emulator] Attaching to existing shared memory segment")
            try:
                self.shm = posix_ipc.SharedMemory(SHM_NAME)
                created_new = False
            except posix_ipc.PermissionsError as e:
                print(f"[ERROR] Cannot access existing shared memory: {e}")
                print(f"[ERROR] Try: ipcrm -m <shmid> to remove orphaned segments")
                print(f"[ERROR] Or: ipcs -m to list all segments")
                raise
        
        try:
            self.buf = mmap.mmap(self.shm.fd, SHM_TOTAL)
            os.close(self.shm.fd)  # Close fd, mmap holds the reference
        except Exception as e:
            print(f"[ERROR] Failed to mmap shared memory: {e}")
            raise

        # Only initialize if we created a new segment
        if created_new:
            # Initialize ring buffer to zeros
            self.buf.seek(0)
            self.buf.write(b'\x00' * SHM_TOTAL)
            self.buf.flush()

            # Write header
            self._write_header()
        else:
            print(f"[emulator] Reusing existing shared memory, skipping init")
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


class DataReader:
    """Read data from various formats (HDF5, .npy, binary, .mcs)."""

    def __init__(self, filepath):
        self.filepath = filepath
        self.file = None
        self.data = None
        self.time_axis = 0
        self.data_key = None
        
        # Auto-detect format and load
        self._load_file()

    def _load_file(self):
        """Auto-detect and load file format."""
        filepath = self.filepath
        
        # Try HDF5 first
        if filepath.endswith(('.h5', '.hdf5', '.HDF5')):
            self._load_hdf5(filepath)
            return
        
        # Try .npy
        if filepath.endswith('.npy'):
            self._load_npy(filepath)
            return
        
        # Try to detect HDF5 by content
        try:
            self._load_hdf5(filepath)
            return
        except:
            pass
        
        # Fallback: try as binary raw file
        print(f"[emulator] WARNING: Unknown format, trying raw binary", file=sys.stderr)
        self._load_raw(filepath)

    def _load_hdf5(self, filepath):
        """Load HDF5 file (handles Maxwell format)."""
        if not HAS_H5PY:
            raise ImportError("h5py required for HDF5 files: pip install h5py")
        
        print(f"[emulator] Opening HDF5 file: {filepath}")
        self.file = h5py.File(filepath, 'r')
        
        # Print top-level structure
        print(f"[emulator] Top-level HDF5 keys: {list(self.file.keys())}")
        
        # Find actual waveform data (2D array with many samples)
        # Maxwell might store data in: recordings/, data_store/, or raw/
        self.data_key = None
        candidates = []
        
        def find_data_arrays(name, obj):
            if isinstance(obj, h5py.Dataset):
                # Look for 2D arrays with substantial size
                if len(obj.shape) == 2 and obj.shape[0] > 100 and obj.shape[1] > 1:
                    candidates.append((name, obj.shape))
                # Also look for very large 1D arrays (might be raveled waveform)
                elif len(obj.shape) == 1 and obj.shape[0] > 10_000:
                    candidates.append((name, obj.shape))
        
        self.file.visititems(find_data_arrays)
        
        if candidates:
            print(f"[emulator] Found candidate data arrays:")
            for name, shape in candidates:
                print(f"  {name}: shape={shape}")
            # Use the largest one
            self.data_key = max(candidates, key=lambda x: np.prod(x[1]))[0]
        
        # If no large arrays found, check common container paths
        if not self.data_key:
            print(f"[emulator] Searching in known container paths...")
            
            # Check recordings/
            if 'recordings' in self.file:
                print(f"[emulator] Found 'recordings' group, exploring...")
                recordings = self.file['recordings']
                for key in recordings.keys():
                    item = recordings[key]
                    if isinstance(item, h5py.Dataset):
                        print(f"  recordings/{key}: shape={item.shape}")
                        if len(item.shape) >= 2 and item.shape[0] > 100:
                            self.data_key = f"recordings/{key}"
                            break
            
            # Check data_store/data0000/ more carefully
            if not self.data_key and 'data_store' in self.file:
                print(f"[emulator] Checking data_store structure...")
                ds = self.file['data_store']
                for well_key in ds.keys():
                    well = ds[well_key]
                    print(f"  data_store/{well_key} keys: {list(well.keys())}")
                    # Look for 'raw', 'data', 'waveform', etc.
                    for data_key in ['raw', 'data', 'waveform', 'continuous']:
                        if data_key in well and isinstance(well[data_key], h5py.Dataset):
                            dset = well[data_key]
                            print(f"    → Found {data_key}: shape={dset.shape}")
                            if len(dset.shape) >= 2 and dset.shape[0] > 100:
                                self.data_key = f"data_store/{well_key}/{data_key}"
                                break
                    if self.data_key:
                        break
        
        if not self.data_key:
            print(f"[emulator] ERROR: No suitable waveform data found!")
            print(f"[emulator] Full HDF5 structure:")
            def print_all(name, obj):
                indent = len(name.split('/')) * 2
                if isinstance(obj, h5py.Dataset):
                    print(f"  {' ' * indent}{name}: shape={obj.shape} dtype={obj.dtype}")
                elif isinstance(obj, h5py.Group):
                    print(f"  {' ' * indent}{name}/ [group]")
            self.file.visititems(print_all)
            
            raise ValueError("No continuous waveform data found in HDF5 file. "
                           "This file may only contain spike data.")
        
        dset = self.file[self.data_key]
        print(f"[emulator] Using dataset '{self.data_key}' shape={dset.shape} dtype={dset.dtype}")
        
        # Determine shape: assume (time, channels) or (channels, time) or (time,) for 1D
        if len(dset.shape) == 1:
            raise ValueError(f"Cannot stream 1D data. Expected 2D array, got shape {dset.shape}")
        elif len(dset.shape) == 2:
            if dset.shape[0] > dset.shape[1]:
                # First dim larger → likely time
                self.time_axis = 0
                self.n_samples, self.n_channels = dset.shape
            else:
                # Second dim larger → likely time
                self.time_axis = 1
                self.n_channels, self.n_samples = dset.shape
        else:
            raise ValueError(f"Expected 2D dataset, got shape {dset.shape}")
        
        print(f"[emulator] Data: {self.n_samples} samples × {self.n_channels} channels")

    def _load_npy(self, filepath):
        """Load .npy file."""
        print(f"[emulator] Opening .npy file: {filepath}")
        self.data = np.load(filepath)
        
        if self.data.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {self.data.shape}")
        
        # Determine shape
        if self.data.shape[0] > self.data.shape[1]:
            self.time_axis = 0
            self.n_samples, self.n_channels = self.data.shape
        else:
            self.time_axis = 1
            self.n_channels, self.n_samples = self.data.shape
        
        print(f"[emulator] Data: {self.n_samples} samples × {self.n_channels} channels")

    def _load_raw(self, filepath):
        """Load raw binary file (assumes float32, must be (time, channels) or (channels, time))."""
        print(f"[emulator] Opening raw binary file: {filepath}")
        # This is a fallback - user would need to specify dimensions
        raise NotImplementedError("Raw binary format requires dimension specification")

    def read_frame(self, sample_index):
        """Read one frame (all channels) at given sample index."""
        if self.file:
            # HDF5 case
            dset = self.file[self.data_key]
            if self.time_axis == 0:
                frame = dset[sample_index, :].astype(np.float32)
            else:
                frame = dset[:, sample_index].astype(np.float32)
        else:
            # In-memory case (.npy, etc.)
            if self.time_axis == 0:
                frame = self.data[sample_index, :].astype(np.float32)
            else:
                frame = self.data[:, sample_index].astype(np.float32)
        return frame

    def get_duration(self, sample_rate=20_000):
        """Get recording duration in seconds."""
        return self.n_samples / sample_rate

    def close(self):
        """Close file if open."""
        if self.file:
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
        # Open recording file (auto-detect format)
        print(f"[emulator] Opening {args.file}...")
        reader = DataReader(args.file)

        # Determine channels to stream (max 8 due to ring buffer size)
        if args.channels:
            channel_ids = args.channels[:8]
        else:
            channel_ids = list(range(min(reader.n_channels, MAX_CHANNELS)))
        
        n_channels = len(channel_ids)

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

