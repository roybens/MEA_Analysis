#!/usr/bin/env python3
"""
setup.py  —  MaxOne Setup + Stimulation
========================================
Loads an existing .cfg file, connects a stimulation unit to your chosen
electrode, and fires biphasic pulses at a fixed amplitude and ISI.

Run this AFTER MaxLab Live Scope is open (mxwserver is already running).
Changes will be reflected live in Scope.

Usage
-----
python setup.py \
    --cfg ~/configs/my_network.cfg \
    --stim-electrode 4887 \
    --amplitude 100 \
    --pulse-width 4 \
    --isi 2.0

Arguments
---------
--cfg               Path to .cfg file saved from Scope
--stim-electrode    Electrode index to stimulate (read from Scope GUI)
--amplitude         Pulse amplitude in bits (1 bit = 2.9 mV, inverting amp)
                    e.g. 100 bits = 290 mV
--pulse-width       Half-phase duration in samples (1 sample = 50 µs)
                    e.g. 4 = 200 µs half-phase = 400 µs total biphasic pulse
--isi               Inter-stimulus interval in seconds
--n-pulses          Total number of pulses to fire (default: infinite, Ctrl-C to stop)
--wells             Wells to use — 0 for MaxOne (default), 0-5 for MaxTwo
"""

import argparse
import time
import sys

import maxlab as mx
import maxlab.chip
import maxlab.util

# ─────────────────────────────────────────────────────────────────────────────
# DEFAULTS  ← edit these if you prefer not to use command-line args
# ─────────────────────────────────────────────────────────────────────────────
DEFAULT_CFG           = ""          # path to your .cfg file
DEFAULT_STIM_ELECTRODE = 4887       # electrode index from Scope
DEFAULT_AMPLITUDE     = 100         # bits (100 × 2.9 mV = 290 mV)
DEFAULT_PULSE_WIDTH   = 4           # samples × 50 µs = 200 µs half-phase
DEFAULT_ISI           = 2.0         # seconds
DEFAULT_N_PULSES      = 0           # 0 = infinite
DEFAULT_WELLS         = [0]         # [0] for MaxOne
# ─────────────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="MaxOne setup + stimulation")
    p.add_argument("--cfg",            type=str,   default=DEFAULT_CFG)
    p.add_argument("--stim-electrode", type=int,   default=DEFAULT_STIM_ELECTRODE)
    p.add_argument("--amplitude",      type=int,   default=DEFAULT_AMPLITUDE)
    p.add_argument("--pulse-width",    type=int,   default=DEFAULT_PULSE_WIDTH)
    p.add_argument("--isi",            type=float, default=DEFAULT_ISI)
    p.add_argument("--n-pulses",       type=int,   default=DEFAULT_N_PULSES,
                   help="0 = run until Ctrl-C")
    p.add_argument("--wells",          type=int,   nargs="+", default=DEFAULT_WELLS)
    return p.parse_args()


def make_biphasic_pulse(amplitude: int, phase_samples: int, dac: int = 0) -> mx.Sequence:
    """
    Single symmetric biphasic pulse.

    DAC range: 0–1023 bits. 512 = 0 V. Inverting amplifier:
        512 - amplitude  →  positive voltage on electrode (cathodic)
        512 + amplitude  →  negative voltage on electrode (anodic)

    amplitude     : bits  (multiply by 2.9 to get mV)
    phase_samples : samples per half-phase (multiply by 50 to get µs)
    """
    seq = mx.Sequence()
    seq.append(mx.chip.DAC(dac, 512 - amplitude))   # cathodic phase
    seq.append(mx.system.DelaySamples(phase_samples))
    seq.append(mx.chip.DAC(dac, 512 + amplitude))   # anodic phase
    seq.append(mx.system.DelaySamples(phase_samples))
    seq.append(mx.chip.DAC(dac, 512))               # return to baseline
    return seq


def setup(cfg_path: str, stim_electrode: int, wells: list) -> int:
    """
    Initialize system, load config, connect stim unit.
    Returns the stim unit index.

    NOTE: We call mx.initialize() which resets the chip to a defined state.
    This will briefly interrupt anything Scope is showing — normal behaviour.
    After this, Scope will reflect the new configuration.
    """
    print("[setup] Initializing system ...")
    mx.util.initialize()
    mx.send(mx.chip.Core().enable_stimulation_power(True))

    # Load routing from .cfg file
    print(f"[setup] Loading config: {cfg_path}")
    array = mx.chip.Array("stimulation")
    array.load_config(cfg_path)

    # Connect stim electrode to a stimulation unit
    array.connect_electrode_to_stimulation(stim_electrode)
    stim_str = array.query_stimulation_at_electrode(stim_electrode)
    if not stim_str:
        print(f"[ERROR] Cannot connect electrode {stim_electrode} to any StimulationUnit.")
        print("        Try the electrode directly adjacent to it in Scope.")
        sys.exit(1)
    unit = int(stim_str)
    print(f"[setup] Electrode {stim_electrode} → StimulationUnit {unit}")

    # Download config to chip and run offset compensation
    array.download(wells)
    mx.util.offset()
    print("[setup] Array downloaded + offset done.")

    # Power up the stim unit in voltage mode
    mx.send(
        mx.chip.StimulationUnit(unit)
          .power_up(True)
          .connect(True)
          .set_voltage_mode()
          .dac_source(0)
    )
    print(f"[setup] StimulationUnit {unit} powered up.")
    return unit


def main():
    args = parse_args()

    if not args.cfg:
        print("[ERROR] --cfg is required. Provide the path to your .cfg file.")
        print("        e.g.  python setup.py --cfg ~/configs/my_network.cfg")
        sys.exit(1)

    amp_mv = args.amplitude * 2.9
    pw_us  = args.pulse_width * 50

    print("=" * 55)
    print("  MaxOne Setup + Stimulation")
    print("=" * 55)
    print(f"  Config file    : {args.cfg}")
    print(f"  Stim electrode : {args.stim_electrode}")
    print(f"  Amplitude      : {args.amplitude} bits  ({amp_mv:.0f} mV)")
    print(f"  Pulse width    : {args.pulse_width} samples  ({pw_us} µs half-phase)")
    print(f"  ISI            : {args.isi} s")
    print(f"  Pulses         : {'infinite (Ctrl-C to stop)' if args.n_pulses == 0 else args.n_pulses}")
    print("=" * 55)

    # ── Setup ─────────────────────────────────────────────────────────────
    stim_unit = setup(args.cfg, args.stim_electrode, args.wells)

    # ── Pulse loop ────────────────────────────────────────────────────────
    print(f"\n[stim] Starting pulses every {args.isi} s  (Ctrl-C to stop)\n")
    pulse_count = 0

    try:
        while True:
            seq = make_biphasic_pulse(args.amplitude, args.pulse_width)
            seq.send()
            pulse_count += 1
            print(f"  pulse {pulse_count}"
                  + (f"/{args.n_pulses}" if args.n_pulses > 0 else "")
                  + f"  amp={args.amplitude} bits ({amp_mv:.0f} mV)")

            if args.n_pulses > 0 and pulse_count >= args.n_pulses:
                break

            time.sleep(args.isi)

    except KeyboardInterrupt:
        pass

    print(f"\n[stim] Done. {pulse_count} pulses fired.")


if __name__ == "__main__":
    main()
