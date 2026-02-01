#!/usr/bin/env python3
"""
checkpoint_dashboard_full.py

Robust Streamlit dashboard for MEA pipeline checkpoint JSON files.

- Fully aligned with ProcessingStage enum (0–8)
- Defensive against missing / partial fields
- Never crashes on malformed checkpoints
- Correct stage, failure, and completeness logic

Run:
  streamlit run checkpoint_dashboard_full.py -- \
    --checkpoint-dir /path/to/checkpoints
"""

import streamlit as st
import json
from pathlib import Path
import pandas as pd
import argparse
from datetime import datetime

# ============================================================
# PIPELINE STAGES (AUTHORITATIVE)
# ============================================================

STAGE_MAP = {
    0: "NOT_STARTED",
    1: "PREPROCESSING",
    2: "PREPROCESSING_COMPLETE",
    3: "SORTING",
    4: "SORTING_COMPLETE",
    5: "ANALYZER",
    6: "ANALYZER_COMPLETE",
    7: "REPORTS",
    8: "REPORTS_COMPLETE",
}

IN_PROGRESS_STAGES = {
    "PREPROCESSING",
    "SORTING",
    "ANALYZER",
    "REPORTS",
}

COMPLETE_STAGES = {
    "PREPROCESSING_COMPLETE",
    "SORTING_COMPLETE",
    "ANALYZER_COMPLETE",
    "REPORTS_COMPLETE",
}

# ============================================================
# SAFE ACCESS
# ============================================================

def safe_get(d, key, default=None):
    try:
        return d.get(key, default)
    except Exception:
        return default

# ============================================================
# ARG PARSING (STREAMLIT-SAFE)
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--checkpoint-dir",
        type=str,
        required=False,
        help="Directory containing checkpoint JSON files",
    )
    known, _ = parser.parse_known_args()
    return known

# ============================================================
# LOAD CHECKPOINTS
# ============================================================

@st.cache_data(ttl=300)
def load_checkpoints_dataframe(checkpoint_dir: str):
    rows = []
    errors = []

    checkpoint_path = Path(checkpoint_dir)
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint dir does not exist: {checkpoint_path}")

    for f in sorted(checkpoint_path.glob("*.json")):
        try:
            raw = json.loads(f.read_text())
        except Exception as e:
            errors.append((str(f), str(e)))
            continue

        # -------------------------
        # Canonical + defensive extraction
        # -------------------------
        project = safe_get(raw, "project_name") or safe_get(raw, "project")
        date = safe_get(raw, "date")
        chip = safe_get(raw, "chip_id") or safe_get(raw, "chip")
        run = safe_get(raw, "run_id")
        well = safe_get(raw, "well_id") or safe_get(raw, "well")
        rec = safe_get(raw, "rec_name")

        stage_num = safe_get(raw, "stage")
        stage_name = STAGE_MAP.get(stage_num, "UNKNOWN")

        failed = (
            bool(safe_get(raw, "failed", False))
            or safe_get(raw, "failed_stage") is not None
            or safe_get(raw, "error") is not None
        )

        failed_stage = safe_get(raw, "failed_stage")
        if failed_stage is not None:
            failed_stage_name = STAGE_MAP.get(failed_stage, failed_stage)
            stage_name = f"FAILED_AT_{failed_stage_name}"

        num_units = (
            safe_get(raw, "num_units_filtered")
            or safe_get(raw, "num_units")
            or safe_get(raw, "n_units")
        )

        analyzer_folder = (
            safe_get(raw, "analyzer_folder")
            or safe_get(raw, "output_dir")
        )

        last_updated = safe_get(raw, "last_updated")

        rows.append({
            "file": f.name,
            "path": str(f.resolve()),
            "project": project,
            "date": date,
            "chip": chip,
            "run": run,
            "well": well,
            "rec": rec,
            "stage": stage_name,
            "stage_num": stage_num,
            "failed": failed,
            "failed_stage": failed_stage,
            "error": safe_get(raw, "error"),
            "num_units": num_units,
            "analyzer_folder": analyzer_folder,
            "last_updated": last_updated,
            "_raw": raw,
        })

    df = pd.DataFrame(rows)

    # -------------------------
    # Normalize missing fields
    # -------------------------
    for col in ["project", "chip", "run", "well", "rec", "stage"]:
        if col in df.columns:
            df[col] = df[col].fillna("—")

    if "num_units" in df.columns:
        df["num_units"] = pd.to_numeric(df["num_units"], errors="coerce")

    if "last_updated" in df.columns:
        df["last_updated"] = pd.to_datetime(df["last_updated"], errors="coerce")

    return df, errors

# ============================================================
# UI HELPERS
# ============================================================

def stage_badge(stage):
    base = stage.replace("FAILED_AT_", "")
    color_map = {
        # completed
        "REPORTS_COMPLETE": "#1f9d55",
        "ANALYZER_COMPLETE": "#1f9d55",
        "SORTING_COMPLETE": "#1f9d55",
        "PREPROCESSING_COMPLETE": "#1f9d55",

        # in progress
        "REPORTS": "#ff9f1c",
        "ANALYZER": "#ff9f1c",
        "SORTING": "#ff9f1c",
        "PREPROCESSING": "#ff9f1c",

        # misc
        "NOT_STARTED": "#95a5a6",
        "UNKNOWN": "#7f8c8d",
    }

    color = color_map.get(base, "#d90429" if stage.startswith("FAILED") else "#2d2d2d")
    label = stage.replace("_", " ")

    return f"""
    <span style="padding:5px 10px;border-radius:6px;
                 background:{color};color:white;font-weight:600;">
        {label}
    </span>
    """

def badge(text, bg):
    return f"<span style='padding:3px 8px;border-radius:6px;background:{bg}'>{text}</span>"

# ============================================================
# STREAMLIT APP
# ============================================================

def run_app(checkpoint_dir):
    st.set_page_config(layout="wide", page_title="MEA Checkpoint Dashboard")
    st.title("📊 MEA Pipeline Checkpoint Dashboard")

    df, errors = load_checkpoints_dataframe(checkpoint_dir)

    # -------------------------
    # Summary tiles
    # -------------------------
    total = len(df)
    completed = (df["stage"] == "REPORTS_COMPLETE").sum()
    failed = df["failed"].sum()
    in_progress = df["stage"].isin(IN_PROGRESS_STAGES).sum()

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total", total)
    c2.metric("Completed", completed)
    c3.metric("In progress", in_progress)
    c4.metric("Failed", failed)

    # -------------------------
    # Filters
    # -------------------------
    st.sidebar.header("Filters")

    for col in ["project", "chip", "run", "well", "stage"]:
        if col not in df.columns:
            continue
        options = sorted(df[col].dropna().unique())
        sel = st.sidebar.multiselect(col, options)
        if sel:
            df = df[df[col].isin(sel)]

    if st.sidebar.checkbox("Failures only"):
        df = df[df["failed"]]

    # -------------------------
    # Table
    # -------------------------
    st.subheader("Checkpoints")

    html = ["<table width='100%' style='border-collapse:collapse'>"]
    html.append("<tr>")
    for c in ["file", "project", "chip", "run", "well", "stage", "num_units", "last_updated"]:
        html.append(f"<th align='left' style='border-bottom:1px solid #ddd'>{c}</th>")
    html.append("</tr>")

    for _, r in df.iterrows():
        ts = r["last_updated"]
        ts_str = ts.strftime("%Y-%m-%d %H:%M:%S") if pd.notna(ts) else "—"

        units_html = badge(int(r["num_units"]), "#f0ad4e") if pd.notna(r["num_units"]) else ""

        html.append("<tr>")
        html.append(f"<td>{r['file']}</td>")
        html.append(f"<td>{r['project']}</td>")
        html.append(f"<td>{r['chip']}</td>")
        html.append(f"<td>{r['run']}</td>")
        html.append(f"<td>{r['well']}</td>")
        html.append(f"<td>{stage_badge(r['stage'])}</td>")
        html.append(f"<td>{units_html}</td>")
        html.append(f"<td>{ts_str}</td>")
        html.append("</tr>")

    html.append("</table>")
    st.markdown("".join(html), unsafe_allow_html=True)

    # -------------------------
    # Inspector
    # -------------------------
    st.subheader("Inspect checkpoint")

    if len(df) == 0:
        st.info("No checkpoints to inspect.")
        return

    sel_file = st.selectbox("Select file", df["file"].tolist())
    row = df[df["file"] == sel_file].iloc[0]

    if row["analyzer_folder"] not in (None, "—"):
        st.text_input("Analyzer / Output folder", row["analyzer_folder"])

    st.json(row["_raw"])

    # -------------------------
    # Parse errors
    # -------------------------
    if errors:
        st.error("Some checkpoint files failed to parse")
        for f, e in errors:
            st.write(f"{f}: {e}")

# ============================================================
# ENTRYPOINT
# ============================================================

if __name__ == "__main__":
    args = parse_args()
    checkpoint_dir = args.checkpoint_dir or "./AnalyzedData/checkpoints"
    run_app(checkpoint_dir)