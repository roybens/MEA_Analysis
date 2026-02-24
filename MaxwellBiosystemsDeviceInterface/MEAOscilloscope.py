#!/usr/bin/env python3
"""
oscilloscope.py  –  MaxOne Live Oscilloscope + Stimulation Controller
======================================================================
Live rolling display of raw MaxOne data via shared memory (mxw_streamer).
Fires stimulation pulses automatically at a fixed ISI.

Controls
--------
  SPACE           freeze / unfreeze the display
  LEFT / RIGHT    scroll backward / forward 1 s through frozen history
  SHIFT+LEFT/RIGHT  fine scroll 100 ms steps
  R               jump back to live (unfreeze + go to present)
  S               save screenshot as PNG
  Q / Escape      quit

Run order
---------
  Terminal 1:  python oscilloscope.py --stim-electrode 4887 \\
                   --detect-electrodes 4886 4885 4888 \\
                   --amplitudes 50 100 150 200 --isi 2.0
  Terminal 2:  ./mxw_streamer 42 55 67 88   (channel indices from MaxLab Live)

Requirements:  pip install numpy matplotlib posix_ipc
"""

import argparse
import mmap
import os
import struct
import threading
import time
import datetime

import numpy as np
import matplotlib
matplotlib.use("TkAgg")        # change to Qt5Agg if TkAgg unavailable
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import maxlab as mx
import maxlab.saving
import maxlab.system

# ─────────────────────────────────────────────────────────────────────────────
# USER PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_STIM_ELECTRODE     = 4887
DEFAULT_DETECT_ELECTRODES  = [4886, 4885, 4888]
DEFAULT_AMPLITUDES         = [50, 100, 150, 200]   # bits
DEFAULT_PULSE_WIDTH        = 4                      # samples × 50 µs
DEFAULT_ISI                = 2.0                   # seconds between pulses
DEFAULT_PULSES_PER_AMP     = 10
DEFAULT_DATA_PATH          = "/tmp"
DEFAULT_TAG                = "scope_session"

# Display
DISPLAY_WINDOW_SEC   = 4.0
DISPLAY_REFRESH_MS   = 50        # 20 fps
SCROLL_STEP_SEC      = 1.0       # coarse (arrow keys)
SCROLL_STEP_FINE_SEC = 0.1       # fine (shift + arrow)
SAMPLING_RATE        = 20_000    # Hz

# Shared memory — must match mxw_streamer.cpp exactly
SHM_NAME     = "/mxw_scope"
SHM_MAGIC    = 0xC0FFEE42
MAX_CHANNELS = 8
RING_FRAMES  = 40_000            # 2 s @ 20 kHz
DATA_OFFSET  = 4096
DATA_SIZE    = RING_FRAMES * MAX_CHANNELS * 4   # float32
STIM_OFFSET  = DATA_OFFSET + DATA_SIZE
SHM_TOTAL    = STIM_OFFSET + RING_FRAMES

HEADER_FMT  = "=II8HQQIF"
HEADER_SIZE = struct.calcsize(HEADER_FMT)

# ─────────────────────────────────────────────────────────────────────────────


def read_header(shm_buf):
    fields = struct.unpack_from(HEADER_FMT, shm_buf, 0)
    return {
        "magic"        : fields[0],
        "n_channels"   : fields[1],
        "channel_ids"  : list(fields[2:10]),
        "write_pos"    : fields[10],
        "frame_number" : fields[11],
        "sample_rate"  : fields[12],
        "lsb_uv"       : fields[13],
    }


def read_ring_slice(shm_buf, end_pos, n_samples, n_channels):
    """
    Return ndarray (n_samples, n_channels) float32 ending at end_pos in the ring.
    """
    end_slot   = int(end_pos) % RING_FRAMES
    start_slot = int(end_pos - n_samples) % RING_FRAMES

    data_view = np.frombuffer(shm_buf, dtype=np.float32,
                               count=RING_FRAMES * MAX_CHANNELS,
                               offset=DATA_OFFSET)
    data_view = data_view.reshape(RING_FRAMES, MAX_CHANNELS)[:, :n_channels]

    if start_slot < end_slot:
        return data_view[start_slot:end_slot].copy()
    else:
        return np.concatenate([data_view[start_slot:],
                                data_view[:end_slot]], axis=0)


# ─────────────────────────────────────────────────────────────────────────────
# MaxOne setup
# ─────────────────────────────────────────────────────────────────────────────

def setup_maxone(stim_el, all_electrodes):
    print("[setup] Initializing MaxOne ...")
    mx.initialize()
    if mx.send(mx.Core().enable_stimulation_power(True)) != "Ok":
        raise RuntimeError("MaxOne init failed.")

    array = mx.Array("stimulation")
    array.reset()
    array.clear_selected_electrodes()
    array.select_electrodes(list(set(all_electrodes + [stim_el])))
    array.select_stimulation_electrodes([stim_el])
    array.route()

    array.connect_electrode_to_stimulation(stim_el)
    stim_str = array.query_stimulation_at_electrode(stim_el)
    if not stim_str:
        raise RuntimeError(f"Cannot connect electrode {stim_el} to any StimulationUnit.")
    unit = int(stim_str)
    print(f"[setup] Electrode {stim_el} → StimulationUnit {unit}")

    array.download()
    mx.util.offset()
    print("[setup] Offset done.")

    mx.send(mx.StimulationUnit(unit)
              .power_up(True).connect(True).set_voltage_mode().dac_source(0))
    print(f"[setup] StimulationUnit {unit} powered up.")
    return array, unit


def make_pulse(amplitude, phase_samples, dac=0):
    seq = mx.Sequence()
    seq.append(mx.DAC(dac, 512 - amplitude))
    seq.append(mx.DelaySamples(phase_samples))
    seq.append(mx.DAC(dac, 512 + amplitude))
    seq.append(mx.DelaySamples(phase_samples))
    seq.append(mx.DAC(dac, 512))
    return seq


# ─────────────────────────────────────────────────────────────────────────────
# Stimulation thread
# ─────────────────────────────────────────────────────────────────────────────

class StimController:
    def __init__(self, amplitudes, pulse_width, isi, pulses_per_amp):
        self.amplitudes     = amplitudes
        self.pulse_width    = pulse_width
        self.isi            = isi
        self.pulses_per_amp = pulses_per_amp
        self.running        = False
        self.current_amp    = amplitudes[0]
        self.amp_idx        = 0
        self.pulse_count    = 0
        self.total_pulses   = 0
        self.stim_log       = []   # list of (wall_time, amp_bits)
        self._lock          = threading.Lock()

    def start(self):
        self.running = True
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def stop(self):
        self.running = False

    def _loop(self):
        while self.running:
            amp = self.amplitudes[self.amp_idx]
            seq = make_pulse(amp, self.pulse_width)
            seq.send()

            with self._lock:
                self.current_amp   = amp
                self.total_pulses += 1
                self.pulse_count  += 1
                self.stim_log.append((time.time(), amp))

            print(f"[stim] amp={amp} bits ({amp*2.9:.0f} mV)  "
                  f"pulse {self.pulse_count}/{self.pulses_per_amp}  "
                  f"total={self.total_pulses}")

            if self.pulse_count >= self.pulses_per_amp:
                self.amp_idx     = (self.amp_idx + 1) % len(self.amplitudes)
                self.pulse_count = 0
                print(f"[stim] → Next amp: {self.amplitudes[self.amp_idx]} bits")

            time.sleep(self.isi)


# ─────────────────────────────────────────────────────────────────────────────
# Oscilloscope
# ─────────────────────────────────────────────────────────────────────────────

class Oscilloscope:
    def __init__(self, n_channels, channel_labels, stim_ctrl, screenshot_dir="/tmp"):
        self.n_ch            = n_channels
        self.labels          = channel_labels
        self.stim_ctrl       = stim_ctrl
        self.screenshot_dir  = screenshot_dir
        self.shm_buf         = None
        self._ready          = False

        self.n_disp = int(DISPLAY_WINDOW_SEC * SAMPLING_RATE)
        self.t_axis = np.linspace(-DISPLAY_WINDOW_SEC, 0, self.n_disp)

        # Freeze / scroll state
        self.frozen           = False
        self.frozen_write_pos = 0
        # scroll_offset: frames back from frozen_write_pos (0 = at freeze moment)
        self.scroll_offset    = 0
        self.max_scroll       = RING_FRAMES - self.n_disp

        self._build_figure()

    def _build_figure(self):
        plt.style.use("dark_background")
        self.fig, self.axes = plt.subplots(
            self.n_ch, 1,
            figsize=(14, 2.4 * self.n_ch),
            sharex=True
        )
        if self.n_ch == 1:
            self.axes = [self.axes]

        colors = ["#00FF41", "#FF6B35", "#FFD700", "#00D4FF",
                  "#FF69B4", "#7FFF00", "#FF4500", "#1E90FF"]
        self.lines       = []
        self.stim_vlines = [[] for _ in range(self.n_ch)]

        for i, ax in enumerate(self.axes):
            line, = ax.plot(self.t_axis, np.zeros(self.n_disp),
                            color=colors[i % len(colors)], lw=0.7,
                            antialiased=False)
            self.lines.append(line)
            ax.set_ylabel(self.labels[i], fontsize=8, color="white")
            ax.set_ylim(-500, 500)
            ax.axhline(0, color="grey", lw=0.4, ls="--", alpha=0.4)
            ax.tick_params(colors="grey", labelsize=7)
            for sp in ax.spines.values():
                sp.set_edgecolor("#333333")
            ax.set_facecolor("#0A0A0A")

        self.axes[-1].set_xlabel("Time (s)", color="white", fontsize=9)
        self.fig.patch.set_facecolor("#050505")

        self.status_text = self.fig.text(
            0.5, 0.985, "Waiting for mxw_streamer...",
            ha="center", va="top", fontsize=9, color="#AAAAAA",
            transform=self.fig.transFigure
        )
        self.freeze_text = self.fig.text(
            0.5, 0.012, "▶ LIVE",
            ha="center", va="bottom", fontsize=9, color="#00FF41",
            transform=self.fig.transFigure, fontweight="bold"
        )
        self.fig.text(
            0.5, 0.001,
            "SPACE freeze/unfreeze  |  ←→ scroll 1s  |  SHIFT+←→ 100ms  |  R→live  |  S screenshot  |  Q quit",
            ha="center", va="bottom", fontsize=6.5, color="#444444",
            transform=self.fig.transFigure
        )

        plt.tight_layout(rect=[0, 0.025, 1, 0.975])
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # ── Key handler ────────────────────────────────────────────────────────

    def _on_key(self, event):
        key = event.key

        if key == " ":
            if not self.frozen:
                if self._ready and self.shm_buf:
                    hdr = read_header(self.shm_buf)
                    self.frozen_write_pos = hdr["write_pos"]
                self.scroll_offset = 0
                self.frozen = True
                print("[scope] ❚❚ FROZEN")
            else:
                self.frozen = False
                self.scroll_offset = 0
                print("[scope] ▶ LIVE")

        elif key == "r":
            self.frozen = False
            self.scroll_offset = 0
            print("[scope] ▶ Back to LIVE")

        elif key in ("left", "shift+left"):
            if self.frozen:
                step = (int(SCROLL_STEP_FINE_SEC * SAMPLING_RATE)
                        if "shift" in key else
                        int(SCROLL_STEP_SEC * SAMPLING_RATE))
                self.scroll_offset = min(self.scroll_offset + step, self.max_scroll)
                print(f"[scope] ← {self.scroll_offset / SAMPLING_RATE:.2f} s before freeze")
            else:
                print("[scope] Freeze first (SPACE), then scroll.")

        elif key in ("right", "shift+right"):
            if self.frozen:
                step = (int(SCROLL_STEP_FINE_SEC * SAMPLING_RATE)
                        if "shift" in key else
                        int(SCROLL_STEP_SEC * SAMPLING_RATE))
                self.scroll_offset = max(self.scroll_offset - step, 0)
                print(f"[scope] → {self.scroll_offset / SAMPLING_RATE:.2f} s before freeze")
            else:
                print("[scope] Freeze first (SPACE), then scroll.")

        elif key == "s":
            ts    = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            fname = os.path.join(self.screenshot_dir, f"scope_{ts}.png")
            self.fig.savefig(fname, dpi=150, bbox_inches="tight",
                             facecolor=self.fig.get_facecolor())
            print(f"[scope] Screenshot → {fname}")

        elif key in ("q", "escape"):
            plt.close("all")

    # ── Shared memory open ─────────────────────────────────────────────────

    def _open_shm(self):
        try:
            import posix_ipc
            shm = posix_ipc.SharedMemory(SHM_NAME)
            self.shm_buf = mmap.mmap(shm.fd, SHM_TOTAL,
                                      mmap.MAP_SHARED, mmap.PROT_READ)
            os.close(shm.fd)
            hdr = read_header(self.shm_buf)
            if hdr["magic"] != SHM_MAGIC:
                self.shm_buf.close()
                self.shm_buf = None
                return False
            self._ready = True
            n = hdr["n_channels"]
            print(f"[scope] Connected — {n} ch: {hdr['channel_ids'][:n]}")
            return True
        except Exception:
            return False

    # ── Stim markers ───────────────────────────────────────────────────────

    def _stim_events_in_view(self, view_end_wall):
        view_start = view_end_wall - DISPLAY_WINDOW_SEC
        with self.stim_ctrl._lock:
            log = list(self.stim_ctrl.stim_log)
        return [(t - view_end_wall, a)
                for t, a in log
                if view_start <= t <= view_end_wall]

    # ── Animation update ───────────────────────────────────────────────────

    def update(self, _frame_num):
        if not self._ready:
            self._open_shm()
            return self.lines + [self.status_text, self.freeze_text]

        hdr = read_header(self.shm_buf)
        if hdr["magic"] != SHM_MAGIC:
            return self.lines + [self.status_text, self.freeze_text]

        live_pos  = hdr["write_pos"]
        n_ch_shm  = min(int(hdr["n_channels"]), self.n_ch)
        lsb_uv    = float(hdr["lsb_uv"])

        if live_pos < self.n_disp:
            return self.lines + [self.status_text, self.freeze_text]

        # Which position to read from
        view_end_pos = (self.frozen_write_pos - self.scroll_offset
                        if self.frozen else live_pos)

        if view_end_pos < self.n_disp:
            return self.lines + [self.status_text, self.freeze_text]

        # Grab data and convert to µV
        raw       = read_ring_slice(self.shm_buf, view_end_pos, self.n_disp, n_ch_shm)
        traces_uv = (raw - 512.0) * lsb_uv

        for i in range(min(self.n_ch, n_ch_shm)):
            trace = traces_uv[:, i]
            self.lines[i].set_ydata(trace)
            p_lo, p_hi = np.percentile(trace, [1, 99])
            margin = max(abs(p_hi - p_lo) * 0.3, 50)
            self.axes[i].set_ylim(p_lo - margin, p_hi + margin)

        # Stim markers
        frames_behind = live_pos - view_end_pos
        view_end_wall = time.time() - frames_behind / SAMPLING_RATE
        stim_events   = self._stim_events_in_view(view_end_wall)

        for i, ax in enumerate(self.axes):
            for vl in self.stim_vlines[i]:
                try:
                    vl.remove()
                except Exception:
                    pass
            self.stim_vlines[i] = []
            for t_rel, amp in stim_events:
                vl = ax.axvline(t_rel, color="#FF3333", lw=1.0, ls="--", alpha=0.85)
                self.stim_vlines[i].append(vl)

        # Status + freeze indicator
        with self.stim_ctrl._lock:
            cur_amp = self.stim_ctrl.current_amp
            total_p = self.stim_ctrl.total_pulses

        self.status_text.set_text(
            f"amp: {cur_amp} bits ({cur_amp*2.9:.0f} mV)  |  "
            f"pulses: {total_p}  |  frame: {hdr['frame_number']}"
        )

        if self.frozen:
            offset_sec = self.scroll_offset / SAMPLING_RATE
            self.freeze_text.set_text(
                f"❚❚ FROZEN  —  {offset_sec:.2f} s before freeze  "
                f"(SPACE=unfreeze  ←→=scroll  R=live)"
            )
            self.freeze_text.set_color("#FF4444")
        else:
            self.freeze_text.set_text("▶ LIVE")
            self.freeze_text.set_color("#00FF41")

        return self.lines + [self.status_text, self.freeze_text]

    def run(self):
        self.ani = animation.FuncAnimation(
            self.fig, self.update,
            interval=DISPLAY_REFRESH_MS,
            blit=False,
            cache_frame_data=False
        )
        plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# CLI + main
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="MaxOne live oscilloscope")
    p.add_argument("--stim-electrode",    type=int,   default=DEFAULT_STIM_ELECTRODE)
    p.add_argument("--detect-electrodes", type=int,   nargs="+",
                   default=DEFAULT_DETECT_ELECTRODES)
    p.add_argument("--amplitudes",        type=int,   nargs="+",
                   default=DEFAULT_AMPLITUDES)
    p.add_argument("--pulse-width",       type=int,   default=DEFAULT_PULSE_WIDTH)
    p.add_argument("--isi",               type=float, default=DEFAULT_ISI)
    p.add_argument("--pulses-per-amp",    type=int,   default=DEFAULT_PULSES_PER_AMP)
    p.add_argument("--data-path",         type=str,   default=DEFAULT_DATA_PATH)
    p.add_argument("--tag",               type=str,   default=DEFAULT_TAG)
    p.add_argument("--no-record",         action="store_true")
    p.add_argument("--screenshot-dir",    type=str,   default="/tmp")
    return p.parse_args()


def main():
    args     = parse_args()
    stim_el  = args.stim_electrode
    det_els  = args.detect_electrodes
    all_els  = list(set([stim_el] + det_els))

    ordered_els = [stim_el] + [e for e in det_els if e != stim_el]
    n_ch   = len(ordered_els)
    labels = [f"el {e}" + (" [stim]" if e == stim_el else " [rec]")
              for e in ordered_els]

    print("=" * 60)
    print("  MaxOne Live Oscilloscope")
    print(f"  Stim electrode   : {stim_el}")
    print(f"  Detect electrodes: {det_els}")
    print(f"  Amplitudes       : {args.amplitudes} bits")
    print(f"  Pulse half-width : {args.pulse_width} × 50µs = {args.pulse_width*50} µs")
    print(f"  ISI              : {args.isi} s  |  {args.pulses_per_amp} pulses/amp")
    print("=" * 60)

    array, stim_unit = setup_maxone(stim_el, all_els)

    saving = None
    if not args.no_record:
        ts       = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        rec_name = f"{args.tag}_{ts}"
        saving   = mx.saving.Saving()
        saving.open_directory(args.data_path)
        saving.group_delete_all()
        saving.group_define(0, "routed")
        saving.start_file(rec_name)
        saving.start_recording([0])
        print(f"[rec]   Recording → {args.data_path}/{rec_name}")

    stim_ctrl = StimController(
        amplitudes     = args.amplitudes,
        pulse_width    = args.pulse_width,
        isi            = args.isi,
        pulses_per_amp = args.pulses_per_amp,
    )

    print("\n[scope] Window opening — start mxw_streamer now, pulses begin in 3 s\n")
    time.sleep(3)
    stim_ctrl.start()

    scope = Oscilloscope(n_ch, labels, stim_ctrl,
                          screenshot_dir=args.screenshot_dir)
    try:
        scope.run()
    finally:
        stim_ctrl.stop()
        if saving:
            saving.stop_recording()
            saving.stop_file()
            saving.group_delete_all()
            print("[rec]   Recording saved.")
        if scope.shm_buf:
            scope.shm_buf.close()
        print("[scope] Done.")


if __name__ == "__main__":
    main()
