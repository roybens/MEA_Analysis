#!/usr/bin/env python3
"""
oscilloscope.py  —  MaxOne Live Oscilloscope Viewer
=====================================================
Reads raw channel data from shared memory written by mxw_streamer (C++)
and displays a live rolling oscilloscope window.

No MaxLab API calls — this script is a pure viewer.
Stimulation is handled by setup.py running separately.

Controls
--------
  SPACE           Freeze / unfreeze
  LEFT / RIGHT    Scroll 1 s steps through frozen history
  SHIFT+LEFT/RIGHT  Fine scroll 100 ms steps
  R               Jump back to live
  S               Save screenshot PNG to --screenshot-dir
  Q / Escape      Quit

Run order
---------
  Terminal 1:  python setup.py --cfg ~/configs/net.cfg --stim-electrode 4887 ...
  Terminal 2:  ./mxw_streamer 42 55 67 88     (channel numbers from Scope GUI)
  Terminal 3:  python oscilloscope.py --channels 42 55 67 88 --labels "stim" "rec1" "rec2" "rec3"

Requirements
------------
  pip install numpy matplotlib posix_ipc
"""

import argparse
import mmap
import os
import struct
import time
import datetime

import numpy as np
import matplotlib
matplotlib.use("TkAgg")        # swap to Qt5Agg if TkAgg unavailable
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ─────────────────────────────────────────────────────────────────────────────
# DISPLAY PARAMETERS  ← adjust to taste
# ─────────────────────────────────────────────────────────────────────────────

DISPLAY_WINDOW_SEC   = 4.0     # seconds of history visible at once
DISPLAY_REFRESH_MS   = 50      # animation interval — 50 ms = 20 fps
SCROLL_STEP_SEC      = 1.0     # coarse scroll (arrow keys)
SCROLL_STEP_FINE_SEC = 0.1     # fine scroll (shift + arrow)
SAMPLING_RATE        = 20_000  # Hz

# ─────────────────────────────────────────────────────────────────────────────
# Shared memory layout — must match mxw_streamer.cpp exactly
# ─────────────────────────────────────────────────────────────────────────────

SHM_NAME     = "/mxw_scope"
SHM_MAGIC    = 0xC0FFEE42
MAX_CHANNELS = 8
RING_FRAMES  = 40_000           # 2 s @ 20 kHz
DATA_OFFSET  = 4096
DATA_SIZE    = RING_FRAMES * MAX_CHANNELS * 4   # float32
SHM_TOTAL    = DATA_OFFSET + DATA_SIZE

# Header format string for struct.unpack (must match ShmHeader in C++)
#   = native byte order, no alignment padding
#   I  magic          uint32
#   I  n_channels     uint32
#   8H channel_ids    uint16 × 8
#   Q  write_pos      uint64
#   Q  frame_number   uint64
#   I  sample_rate    uint32
#   f  lsb_uv         float32
HEADER_FMT  = "=II8HQQIF"
HEADER_SIZE = struct.calcsize(HEADER_FMT)

# ─────────────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="MaxOne live oscilloscope viewer")
    p.add_argument("--channels", type=int, nargs="+", required=False,
                   help="Channel numbers to display (from Scope GUI). "
                        "Must match what mxw_streamer was started with.")
    p.add_argument("--labels", type=str, nargs="+", default=None,
                   help="Labels for each channel row, e.g. 'stim' 'rec1' 'rec2'")
    p.add_argument("--window", type=float, default=DISPLAY_WINDOW_SEC,
                   help="Display window in seconds (default 4.0)")
    p.add_argument("--screenshot-dir", type=str, default="/tmp",
                   help="Directory for screenshots saved with S key")
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Shared memory helpers
# ─────────────────────────────────────────────────────────────────────────────

def read_header(shm_buf):
    vals = struct.unpack_from(HEADER_FMT, shm_buf, 0)
    return {
        "magic"       : vals[0],
        "n_channels"  : vals[1],
        "channel_ids" : list(vals[2:10]),
        "write_pos"   : vals[10],
        "frame_number": vals[11],
        "sample_rate" : vals[12],
        "lsb_uv"      : vals[13],
    }


def read_ring_slice(shm_buf, end_pos, n_samples, n_channels):
    """
    Read n_samples frames ending at end_pos from the ring buffer.
    Returns ndarray shape (n_samples, n_channels) in raw ADC float32 values.
    end_pos is absolute (not modded); we mod internally.
    """
    end_slot   = int(end_pos) % RING_FRAMES
    start_slot = int(end_pos - n_samples) % RING_FRAMES

    # View the data region as a (RING_FRAMES × MAX_CHANNELS) float32 array
    arr = np.frombuffer(shm_buf, dtype=np.float32,
                         count=RING_FRAMES * MAX_CHANNELS,
                         offset=DATA_OFFSET
                         ).reshape(RING_FRAMES, MAX_CHANNELS)[:, :n_channels]

    if start_slot < end_slot:
        return arr[start_slot:end_slot].copy()
    else:
        # Wraps around the ring
        return np.concatenate([arr[start_slot:], arr[:end_slot]], axis=0)


# ─────────────────────────────────────────────────────────────────────────────
# Oscilloscope class
# ─────────────────────────────────────────────────────────────────────────────

class Oscilloscope:

    COLORS = ["#00FF41", "#FF6B35", "#FFD700", "#00D4FF",
              "#FF69B4", "#7FFF00", "#FF4500", "#1E90FF"]

    def __init__(self, n_channels, labels, window_sec, screenshot_dir):
        self.n_ch           = n_channels
        self.labels         = labels
        self.window_sec     = window_sec
        self.screenshot_dir = screenshot_dir

        self.n_disp = int(window_sec * SAMPLING_RATE)
        self.t_axis = np.linspace(-window_sec, 0, self.n_disp)

        # Shared memory state
        self.shm_buf  = None
        self._ready   = False

        # Freeze / scroll state
        self.frozen           = False
        self.frozen_write_pos = 0
        self.scroll_offset    = 0          # frames back from frozen_write_pos
        self.max_scroll       = RING_FRAMES - self.n_disp

        self._build_figure()

    # ── Figure ────────────────────────────────────────────────────────────

    def _build_figure(self):
        plt.style.use("dark_background")
        self.fig, axes = plt.subplots(
            self.n_ch, 1,
            figsize=(14, 2.4 * self.n_ch),
            sharex=True
        )
        self.axes = [axes] if self.n_ch == 1 else list(axes)

        self.lines = []
        for i, ax in enumerate(self.axes):
            line, = ax.plot(self.t_axis, np.zeros(self.n_disp),
                            color=self.COLORS[i % len(self.COLORS)],
                            lw=0.7, antialiased=False)
            self.lines.append(line)
            ax.set_ylabel(self.labels[i], fontsize=8, color="white")
            ax.set_ylim(-500, 500)
            ax.axhline(0, color="#333333", lw=0.5, ls="--")
            ax.tick_params(colors="grey", labelsize=7)
            for sp in ax.spines.values():
                sp.set_edgecolor("#222222")
            ax.set_facecolor("#080808")

        self.axes[-1].set_xlabel("Time (s)", color="white", fontsize=9)
        self.fig.patch.set_facecolor("#030303")

        # Status line — top
        self.status_text = self.fig.text(
            0.5, 0.987, "Waiting for mxw_streamer...",
            ha="center", va="top", fontsize=9, color="#888888",
            transform=self.fig.transFigure
        )

        # Live / frozen indicator — bottom centre
        self.state_text = self.fig.text(
            0.5, 0.013, "▶ LIVE",
            ha="center", va="bottom", fontsize=9, color="#00FF41",
            fontweight="bold", transform=self.fig.transFigure
        )

        # Key hint — very bottom
        self.fig.text(
            0.5, 0.001,
            "SPACE freeze  ·  ←→ scroll 1s  ·  SHIFT+←→ 100ms  ·  R live  ·  S screenshot  ·  Q quit",
            ha="center", va="bottom", fontsize=6.5, color="#333333",
            transform=self.fig.transFigure
        )

        plt.tight_layout(rect=[0, 0.025, 1, 0.975])
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # ── Key handler ───────────────────────────────────────────────────────

    def _on_key(self, event):
        k = event.key

        if k == " ":
            if not self.frozen:
                if self._ready:
                    self.frozen_write_pos = read_header(self.shm_buf)["write_pos"]
                self.scroll_offset = 0
                self.frozen = True
                print("[scope] ❚❚ FROZEN")
            else:
                self.frozen = False
                self.scroll_offset = 0
                print("[scope] ▶ LIVE")

        elif k == "r":
            self.frozen = False
            self.scroll_offset = 0
            print("[scope] ▶ Back to LIVE")

        elif k in ("left", "shift+left"):
            if not self.frozen:
                print("[scope] Press SPACE to freeze first, then scroll.")
                return
            step = (int(SCROLL_STEP_FINE_SEC * SAMPLING_RATE) if "shift" in k
                    else int(SCROLL_STEP_SEC * SAMPLING_RATE))
            self.scroll_offset = min(self.scroll_offset + step, self.max_scroll)
            print(f"[scope] ← {self.scroll_offset / SAMPLING_RATE:.2f} s before freeze")

        elif k in ("right", "shift+right"):
            if not self.frozen:
                print("[scope] Press SPACE to freeze first, then scroll.")
                return
            step = (int(SCROLL_STEP_FINE_SEC * SAMPLING_RATE) if "shift" in k
                    else int(SCROLL_STEP_SEC * SAMPLING_RATE))
            self.scroll_offset = max(self.scroll_offset - step, 0)
            print(f"[scope] → {self.scroll_offset / SAMPLING_RATE:.2f} s before freeze")

        elif k == "s":
            ts    = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            fname = os.path.join(self.screenshot_dir, f"scope_{ts}.png")
            self.fig.savefig(fname, dpi=150, bbox_inches="tight",
                             facecolor=self.fig.get_facecolor())
            print(f"[scope] Screenshot saved → {fname}")

        elif k in ("q", "escape"):
            plt.close("all")

    # ── Shared memory connect ─────────────────────────────────────────────

    def _try_connect(self):
        try:
            import posix_ipc
            shm = posix_ipc.SharedMemory(SHM_NAME)
            buf = mmap.mmap(shm.fd, SHM_TOTAL, mmap.MAP_SHARED, mmap.PROT_READ)
            os.close(shm.fd)
            hdr = read_header(buf)
            if hdr["magic"] != SHM_MAGIC:
                buf.close()
                return
            self.shm_buf = buf
            self._ready  = True
            n   = hdr["n_channels"]
            ids = hdr["channel_ids"][:n]
            print(f"[scope] Connected to mxw_streamer — {n} channels: {ids}")
            # If no labels were provided, auto-label from channel IDs
            if len(self.labels) < self.n_ch:
                self.labels = [f"ch {ids[i]}" for i in range(self.n_ch)]
                for i, ax in enumerate(self.axes):
                    ax.set_ylabel(self.labels[i], fontsize=8, color="white")
        except Exception:
            pass   # streamer not started yet — will retry next frame

    # ── Animation update (called every DISPLAY_REFRESH_MS) ───────────────

    def update(self, _):
        if not self._ready:
            self._try_connect()
            return self.lines + [self.status_text, self.state_text]

        hdr = read_header(self.shm_buf)
        if hdr["magic"] != SHM_MAGIC:
            return self.lines + [self.status_text, self.state_text]

        live_pos = hdr["write_pos"]
        n_ch_shm = min(int(hdr["n_channels"]), self.n_ch)
        lsb_uv   = float(hdr["lsb_uv"])

        if live_pos < self.n_disp:
            # Not enough data yet
            return self.lines + [self.status_text, self.state_text]

        # Which position to read from
        view_end = (self.frozen_write_pos - self.scroll_offset
                    if self.frozen else live_pos)

        if view_end < self.n_disp:
            return self.lines + [self.status_text, self.state_text]

        # Pull data from ring buffer and convert to µV
        raw       = read_ring_slice(self.shm_buf, view_end, self.n_disp, n_ch_shm)
        traces_uv = (raw - 512.0) * lsb_uv   # shape: (n_disp, n_ch_shm)

        # Update each channel trace
        for i in range(min(self.n_ch, n_ch_shm)):
            trace = traces_uv[:, i]
            self.lines[i].set_ydata(trace)

            # Auto-scale y axis
            lo, hi = np.percentile(trace, [1, 99])
            margin = max(abs(hi - lo) * 0.3, 50)   # minimum ±50 µV margin
            self.axes[i].set_ylim(lo - margin, hi + margin)

        # Status text
        self.status_text.set_text(
            f"frame {hdr['frame_number']}  ·  "
            f"{hdr['n_channels']} channels  ·  "
            f"{'FROZEN' if self.frozen else 'LIVE'}"
        )

        if self.frozen:
            offset_s = self.scroll_offset / SAMPLING_RATE
            self.state_text.set_text(
                f"❚❚ FROZEN  —  viewing {offset_s:.2f} s before freeze  "
                f"·  SPACE=unfreeze  ·  ←→=scroll  ·  R=live"
            )
            self.state_text.set_color("#FF4444")
        else:
            self.state_text.set_text("▶ LIVE")
            self.state_text.set_color("#00FF41")

        return self.lines + [self.status_text, self.state_text]

    # ── Run ───────────────────────────────────────────────────────────────

    def run(self):
        self.ani = animation.FuncAnimation(
            self.fig, self.update,
            interval=DISPLAY_REFRESH_MS,
            blit=False,
            cache_frame_data=False
        )
        plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # Determine number of channels and labels
    if args.channels:
        n_ch   = len(args.channels)
        labels = args.labels if args.labels else [f"ch {c}" for c in args.channels]
    else:
        # No channels specified — will auto-detect from shared memory on connect
        n_ch   = MAX_CHANNELS
        labels = []

    # Pad or trim labels to match n_ch
    while len(labels) < n_ch:
        labels.append(f"ch {len(labels)}")
    labels = labels[:n_ch]

    window_sec = args.window

    print("=" * 50)
    print("  MaxOne Oscilloscope Viewer")
    print("=" * 50)
    print(f"  Channels    : {args.channels or 'auto-detect from streamer'}")
    print(f"  Labels      : {labels[:n_ch]}")
    print(f"  Window      : {window_sec} s")
    print(f"  Screenshots : {args.screenshot_dir}")
    print("=" * 50)
    print("\n  Make sure mxw_streamer is running.")
    print("  The display will start as soon as it connects.\n")

    scope = Oscilloscope(
        n_channels     = n_ch,
        labels         = labels,
        window_sec     = window_sec,
        screenshot_dir = args.screenshot_dir,
    )
    try:
        scope.run()
    finally:
        if scope.shm_buf:
            scope.shm_buf.close()
        print("[scope] Done.")


if __name__ == "__main__":
    main()
