# MaxOne Oscilloscope System — HOWTO

## Files

| File | Language | Does |
|------|----------|------|
| `setup.py` | Python | Init chip, load .cfg, power stim unit, fire pulses |
| `mxw_streamer.cpp` | C++ | Tap raw stream → shared memory ring buffer |
| `oscilloscope.py` | Python | Read shared memory → live display (no MaxLab API) |
| `recorder.py` | Python | Block record all channels via maxlab.saving API |
| `emulator.py` | Python | Read HDF5 file, emulate streaming to shared memory |

---

## Compile the streamer

```bash
cd ~/MaxLab/share
unzip libmaxlab-*.zip     # skip if already done
cd maxlab_lib

g++ -std=c++17 -O3 -o mxw_streamer /path/to/mxw_streamer.cpp \
    -I$(pwd)/include \
    -L$(pwd)/lib \
    -lmaxlab -lrt -lpthread \
    -Wl,-rpath,$(pwd)/lib
```

---

## Install Python dependencies

```bash
pip install numpy matplotlib posix_ipc h5py
```

For the original C++ streamer approach, `h5py` is optional.

---

## Run order (three terminals)

### Terminal 1 — Setup + stimulation
```bash
python setup.py \
    --cfg ~/configs/my_network.cfg \
    --stim-electrode 4887 \
    --amplitude 100 \
    --pulse-width 4 \
    --isi 2.0
```

| Parameter | Meaning |
|-----------|---------|
| `--cfg` | .cfg file saved from Scope |
| `--stim-electrode` | electrode index — read from Scope GUI |
| `--amplitude` | bits (1 bit = 2.9 mV). Start low (50–100) |
| `--pulse-width` | half-phase samples (× 50 µs). 4 = 200 µs |
| `--isi` | seconds between pulses |
| `--n-pulses` | omit for infinite (Ctrl-C to stop) |

### Terminal 2 — C++ streamer
Read channel numbers from MaxLab Live Scope GUI (hover over electrode → shows channel number).
```bash
./mxw_streamer 42 55 67 88
```
Order matters — first channel = top row in oscilloscope.

### Terminal 3 — Oscilloscope viewer
```bash
python oscilloscope.py \
    --channels 42 55 67 88 \
    --labels "stim el4887" "rec el4886" "rec el4885" "rec el4888"
```
`--channels` and `--labels` are optional — if omitted, auto-detected from the streamer.

---

## Oscilloscope controls

| Key | Action |
|-----|--------|
| `SPACE` | Freeze / unfreeze |
| `←` / `→` | Scroll 1 s backward / forward through frozen history |
| `SHIFT+←` / `SHIFT+→` | Fine scroll 100 ms steps |
| `R` | Jump back to live |
| `S` | Save screenshot PNG |
| `Q` / `Esc` | Quit |

---

## Amplitude tuning workflow

1. Start with `--amplitude 50` (≈145 mV) — you'll see only the artifact
2. Freeze (SPACE) after a pulse, scroll back to inspect the waveform
3. Look on the recording channels for a secondary waveform 2–15 ms after the pulse — that's the evoked AP
4. Ctrl-C setup.py, increase `--amplitude`, re-run setup.py
5. Repeat until AP appears reliably on the detection channel

---

## What to expect on screen

```
Row 0  stim el4887 [green]
  ──────────────⋀⋁──────────────────   ← large fast artifact at stim time

Row 1  rec el4886 [orange]
  ────────────────────⋀────────────   ← evoked AP ~5 ms after stim (smaller, slower)

Row 2  rec el4885 [yellow]
  ──────────────────────────────────   ← nothing if too far / no AP

Row 3  rec el4888 [cyan]
  ──────────────────────────────────
```

Artifact: large amplitude, fast decay, starts exactly at stim pulse.
Evoked AP: ~50–300 µV, slower shape, appears 2–15 ms post-stim on neighbor channels.

---

## Python-Only Alternative (No C++ Compilation Required)

If `mxw_streamer` fails to compile due to library ABI mismatches, use the Python POC approach:
**recorder.py** captures data directly via MaxLab API, then **emulator.py** replays it to shared memory.

### Files

| File | Does |
|------|------|
| `recorder.py` | Block record all channels to HDF5 file via `maxlab.saving` API |
| `emulator.py` | Read HDF5 file, stream to shared memory at real-time speed |

### Install extra Python dependencies

```bash
pip install h5py
```

### Run order (four terminals)

This approach sacrifices real-time display (~30–60 sec latency) for complete Python-only operation.

#### Terminal 1 — Setup + stimulation (unchanged)
```bash
python setup.py \
    --cfg ~/configs/my_network.cfg \
    --stim-electrode 4887 \
    --amplitude 100 \
    --pulse-width 4 \
    --isi 2.0 \
    --n-pulses 10
```
Fires 10 biphasic pulses at your chosen parameters.

#### Terminal 2 — Start block recording
Start this **while Terminal 1 is running** to capture the stimulus:
```bash
python recorder.py \
    --duration 20 \
    --output /tmp/mea_recording.h5 \
    --wells 0
```

Waits 20 seconds and saves all channels to `/tmp/mea_recording.h5`.

#### Terminal 3 — Oscilloscope viewer (unchanged)
```bash
python oscilloscope.py \
    --channels 4887 4886 4885 4888 \
    --labels "stim el4887" "rec el4886" "rec el4885" "rec el4888"
```

Will show "Waiting for mxw_scope..." until emulator starts.

#### Terminal 4 — Start emulator (after recorder finishes)
```bash
python emulator.py /tmp/mea_recording.h5 \
    --speed 1.0
```

Reads HDF5 file and streams to shared memory. Oscilloscope displays live-ish traces (replayed at 20 kHz).

### Workflow for parameter tuning

1. **Initial setup**: Run Terminals 1–3
2. **Record stimulus**: Run Terminal 2 with chosen `--amplitude`, `--pulse-width`, `--isi`
3. **Replay data**: Once Terminal 2 finishes, run Terminal 4
4. **Inspect waveforms**: Use oscilloscope controls (SPACE, ←/→, S) to examine evoked APs
5. **Adjust parameters**: If AP isn't visible, Ctrl-C Terminals 1–4, change `--amplitude`, repeat

### Advantages & Limitations

**Advantages:**
- No C++ compilation needed
- Full Python workflow via maxlab API
- Complete data capture (all 26,400 channels)
- Post-hoc visualization identical to real-time oscilloscope

**Limitations:**
- ~30–60 second delay (download time + playback)
- Cannot adjust amplitude during recording (must restart)
- Requires h5py (`pip install h5py`)

---

## Troubleshooting

**"Waiting for mxw_streamer..."** — streamer not yet started or wrong SHM_NAME.
Start `./mxw_streamer` and the display connects automatically.

**Flat traces** — wrong channel numbers. Double-check in Scope GUI.

**`posix_ipc` import error** — `pip install posix_ipc`

**TkAgg error** — change `matplotlib.use("TkAgg")` to `"Qt5Agg"` in oscilloscope.py
and `pip install PyQt5`

**setup.py: "Cannot connect electrode"** — try the adjacent electrode index.
Routing constraints sometimes prevent a specific electrode from connecting to a stim unit.
