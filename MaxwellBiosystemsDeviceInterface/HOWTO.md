# MaxOne Oscilloscope System — HOWTO

## Files

| File | Language | Does |
|------|----------|------|
| `setup.py` | Python | Init chip, load .cfg, power stim unit, fire pulses |
| `mxw_streamer.cpp` | C++ | Tap raw stream → shared memory ring buffer |
| `oscilloscope.py` | Python | Read shared memory → live display (no MaxLab API) |

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
pip install numpy matplotlib posix_ipc
```

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

## Troubleshooting

**"Waiting for mxw_streamer..."** — streamer not yet started or wrong SHM_NAME.
Start `./mxw_streamer` and the display connects automatically.

**Flat traces** — wrong channel numbers. Double-check in Scope GUI.

**`posix_ipc` import error** — `pip install posix_ipc`

**TkAgg error** — change `matplotlib.use("TkAgg")` to `"Qt5Agg"` in oscilloscope.py
and `pip install PyQt5`

**setup.py: "Cannot connect electrode"** — try the adjacent electrode index.
Routing constraints sometimes prevent a specific electrode from connecting to a stim unit.
