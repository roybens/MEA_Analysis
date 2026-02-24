#!/usr/bin/env python3
"""
recorder.py  —  MaxLab Block Recording via Python API
=======================================================
Records all channels to disk using maxlab.saving.Saving() API.

Run this AFTER setup.py has started stimulation.

Usage
-----
python recorder.py --duration 20 --output /tmp/mea_recording.h5 --wells 0

Arguments
---------
--duration    Duration of recording in seconds (default: 20)
--output      Output file path (default: /tmp/mea_recording.h5)
--wells       Wells to record from, comma-separated (default: 0)
"""

import argparse
import time
import sys
import os
from datetime import datetime

import maxlab.saving
import maxlab.util


def parse_args():
    p = argparse.ArgumentParser(description="MaxLab block recording to disk")
    p.add_argument("--duration", type=int, default=20,
                   help="Recording duration in seconds (default: 20)")
    p.add_argument("--output", type=str, default="/tmp/mea_recording.h5",
                   help="Output file path (default: /tmp/mea_recording.h5)")
    p.add_argument("--wells", type=str, default="0",
                   help="Wells to record (comma-separated, default: 0)")
    return p.parse_args()


def main():
    args = parse_args()

    # Parse wells
    wells = [int(w.strip()) for w in args.wells.split(",")]

    print("=" * 60)
    print("  MaxLab Block Recording")
    print("=" * 60)
    print(f"  Duration    : {args.duration} s")
    print(f"  Output file : {args.output}")
    print(f"  Wells       : {wells}")
    print("=" * 60)
    print()

    try:
        # Ensure output directory exists
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # Get just the filename (no directory)
        filename_only = os.path.splitext(os.path.basename(args.output))[0]

        print(f"[recorder] Initializing maxlab...")
        saving = maxlab.saving.Saving()

        print(f"[recorder] Opening directory: {output_dir if output_dir else '/tmp'}")
        saving.open_directory(output_dir if output_dir else "/tmp")

        print(f"[recorder] Creating file: {filename_only}")
        saving.start_file(filename_only)

        print(f"[recorder] Starting recording on wells {wells} for {args.duration} seconds...")
        saving.start_recording(wells=wells)

        # Record for specified duration
        for i in range(args.duration):
            print(f"  [{i+1}/{args.duration}]", end="\r")
            time.sleep(1)
        print()

        print(f"[recorder] Stopping recording...")
        saving.stop_recording()

        print(f"[recorder] Closing file...")
        saving.stop_file()

        print()
        print(f"[recorder] ✓ Recording complete!")
        print(f"[recorder] Output file: {args.output}")
        print()

    except KeyboardInterrupt:
        print("\n[recorder] Interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
