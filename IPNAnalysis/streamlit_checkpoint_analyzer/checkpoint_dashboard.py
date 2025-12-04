#!/usr/bin/env python3
"""
checkpoint_dashboard_full.py

Improved Streamlit dashboard for MEA pipeline checkpoint JSON files.

Usage:
    streamlit run checkpoint_dashboard_full.py -- --checkpoint-dir /path/to/AnalyzedData/checkpoints

Features:
 - Loads all JSON checkpoints from a directory
 - Stage-aware summary tiles (completed / failed / in-progress / total)
 - Run-level completeness table (completed wells / total wells)
 - Filters: project, date, chip, run, well, stage, failed-only
 - Color-coded stage badges and unit-count badge
 - JSON inspector with extensions expander
 - Export filtered view to CSV
"""

import streamlit as st
import json
from pathlib import Path
import pandas as pd
import argparse
import sys
import textwrap
from datetime import datetime

# -------------------------
# Helpers
# -------------------------
def parse_args():
    """
    Accepts --checkpoint-dir via argv after -- (streamlit passes the rest)
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--checkpoint-dir", type=str, required=False,
                        help="Path to checkpoint directory (default: ./AnalyzedData/checkpoints)")
    known, _ = parser.parse_known_args()
    return known

@st.cache_data(ttl=300)
def load_checkpoints_dataframe(checkpoint_dir: str):
    """
    Loads all .json checkpoint files into a DataFrame.
    Returns (df, load_errors_list)
    """
    checkpoint_path = Path(checkpoint_dir)
    rows = []
    errors = []

    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint dir does not exist: {checkpoint_path}")

    files = sorted(checkpoint_path.glob("*.json"))

    for f in files:
        try:
            raw = json.loads(f.read_text())
        except Exception as e:
            errors.append((str(f), str(e)))
            # Add a minimal row for visibility
            rows.append({
                "file": f.name,
                "path": str(f.resolve()),
                "project": None,
                "date": None,
                "chip": None,
                "run": None,
                "well": None,
                "rec": None,
                "stage": "PARSE_ERROR",
                "stage_num": None,
                "failed": True,
                "error": f"parse_error:{e}",
                "extensions_computed": [],
                "num_units_filtered": None,
                "last_updated": None,
                "_raw": {}
            })
            continue

        # Defensive extraction with defaults
        project = raw.get("project_name") or raw.get("project") or None
        date = raw.get("date")
        chip = raw.get("chip_id") or raw.get("chip")
        run = raw.get("run_id")
        well = raw.get("well_id")
        rec = raw.get("rec_name")
        stage_name = raw.get("stage_name") or raw.get("stage") or None
        # If stage numeric present convert to name mapping if possible
        stage_num = raw.get("stage") if isinstance(raw.get("stage"), (int, float)) else None
        failed = raw.get("failed", False) or (stage_name == "FAILED")
        error = raw.get("error")
        exts = raw.get("extensions_computed") or raw.get("extensions") or []
        num_units = raw.get("num_units_filtered") or raw.get("num_units") or None
        last_updated = raw.get("last_updated")

        # Normalize stage_name string if numeric stage provided
        if stage_name is None and stage_num is not None:
            mapping = {
                -1: "FAILED",
                0: "NOT_STARTED",
                1: "SORTING_COMPLETE",
                2: "ANALYZER_COMPLETE",
                3: "REPORTS_COMPLETE"
            }
            stage_name = mapping.get(int(stage_num), str(stage_num))

        # Fallback stage
        if stage_name is None:
            stage_name = "UNKNOWN"

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
            "failed": bool(failed),
            "error": error,
            "extensions_computed": exts,
            "num_units_filtered": num_units,
            "last_updated": last_updated,
            "_raw": raw
        })

    df = pd.DataFrame(rows)

    # Normalize dtypes and fillna
    if "last_updated" in df.columns:
        def try_parse(dt):
            try:
                return pd.to_datetime(dt)
            except Exception:
                return pd.NaT
        df["last_updated"] = df["last_updated"].apply(try_parse)

    return df, errors

def stage_color_html(stage):
    """
    Return HTML colored badge for a stage.
    """
    color_map = {
        "REPORTS_COMPLETE": "#1f9d55",   # green
        "ANALYZER_COMPLETE": "#ff9f1c",  # orange
        "SORTING_COMPLETE": "#2874a6",   # blue
        "NOT_STARTED": "#95a5a6",        # grey
        "FAILED": "#d90429",             # red
        "PARSE_ERROR": "#d90429",
        "UNKNOWN": "#7f8c8d"
    }
    color = color_map.get(stage, "#2d2d2d")
    display = stage.replace("_", " ").title()
    html = f"""
    <span style="display:inline-block;padding:6px 10px;
                 border-radius:6px;background:{color};color:white;
                 font-weight:600;font-size:0.9em;">
        {display}
    </span>
    """
    return html

def small_badge(text, bg="#e0e0e0"):
    return f"""<span style="display:inline-block;padding:3px 8px;border-radius:6px;background:{bg};font-size:0.85em">{text}</span>"""

def safe_str(x):
    return "" if x is None else str(x)

# -------------------------
# Streamlit UI
# -------------------------
def run_app(checkpoint_dir):
    st.set_page_config(page_title="MEA Checkpoint Dashboard", layout="wide")
    st.title("📊 MEA Pipeline Checkpoint Dashboard — Full")

    # Load data
    try:
        df, errors = load_checkpoints_dataframe(checkpoint_dir)
    except FileNotFoundError as e:
        st.error(str(e))
        st.stop()

    # Top summary tiles
    total = len(df)
    completed = int((df["stage"] == "REPORTS_COMPLETE").sum())
    failed = int(df["failed"].sum())
    in_progress = int(((df["stage"] != "REPORTS_COMPLETE") & (~df["failed"])).sum())
    unknown = int((df["stage"] == "UNKNOWN").sum())

    col1, col2, col3, col4, col5 = st.columns([1,1,1,1,1])
    col1.metric("Total checkpoints", total)
    col2.metric("Completed (reports)", completed)
    col3.metric("Failed", failed)
    col4.metric("In progress", in_progress)
    col5.metric("Unknown", unknown)

    st.markdown("---")

    # Sidebar filters
    st.sidebar.header("Filters / Query")
    st.sidebar.write("Narrow down which checkpoints to show")

    projects = sorted(df["project"].dropna().unique().tolist())
    chips = sorted(df["chip"].dropna().unique().tolist())
    runs = sorted([int(x) for x in df["run"].dropna().unique().tolist()]) if not df["run"].dropna().empty else []
    wells = sorted(df["well"].dropna().unique().tolist())
    stages = sorted(df["stage"].dropna().unique().tolist())

    sel_projects = st.sidebar.multiselect("Project", options=projects)
    sel_chips = st.sidebar.multiselect("Chip ID", options=chips)
    sel_runs = st.sidebar.multiselect("Run ID", options=runs)
    sel_wells = st.sidebar.multiselect("Well", options=wells)
    sel_stages = st.sidebar.multiselect("Stage", options=stages)
    failed_only = st.sidebar.checkbox("Show only failures", value=False)
    min_units = st.sidebar.number_input("Min units filtered (leave 0)", min_value=0, value=0, step=1)

    # Apply filters
    view_df = df.copy()
    if sel_projects:
        view_df = view_df[view_df["project"].isin(sel_projects)]
    if sel_chips:
        view_df = view_df[view_df["chip"].isin(sel_chips)]
    if sel_runs:
        # our runs might be ints or strings; normalize both sides to int if possible
        def to_int_safe(x):
            try:
                return int(x)
            except:
                return None
        sel_runs_int = [int(x) for x in sel_runs]
        view_df = view_df[view_df["run"].apply(lambda x: (to_int_safe(x) in sel_runs_int))]
    if sel_wells:
        view_df = view_df[view_df["well"].isin(sel_wells)]
    if sel_stages:
        view_df = view_df[view_df["stage"].isin(sel_stages)]
    if failed_only:
        view_df = view_df[view_df["failed"] == True]
    if min_units > 0:
        view_df = view_df[view_df["num_units_filtered"].apply(lambda x: (x is not None) and (x >= min_units))]

    st.write(f"Showing **{len(view_df)}** checkpoints (filtered).")

    # Export filtered CSV
    csv_btn_col1, csv_btn_col2 = st.columns([1,3])
    if csv_btn_col1.button("Export CSV"):
        tmp_csv = "checkpoint_filtered_export.csv"
        view_df_out = view_df.drop(columns=["_raw", "extensions_computed"])
        view_df_out.to_csv(tmp_csv, index=False)
        st.success(f"Saved CSV -> {tmp_csv}")
        st.download_button("Download CSV", data=open(tmp_csv, "rb"), file_name=tmp_csv)

    # Run-level completeness
    st.subheader("Run-level completeness")
    if "run" in view_df.columns and not view_df["run"].dropna().empty:
        run_group = view_df.groupby("run").agg(
            wells_total=("well", "count"),
            wells_completed=("stage", lambda s: (s == "REPORTS_COMPLETE").sum())
        )
        run_group["complete_fraction"] = (run_group["wells_completed"] / run_group["wells_total"]).map("{:.2f}".format)
        st.dataframe(run_group.reset_index(), use_container_width=True)
    else:
        st.info("No run-level data available in the filtered selection.")

    st.markdown("---")
    st.subheader("Checkpoint Table")

    # Build display table: customize columns and small badges
    display_cols = ["file", "project", "date", "chip", "run", "well", "rec", "stage", "failed", "num_units_filtered", "last_updated"]
    display_df = view_df[display_cols].copy()
    # Format date
    display_df["last_updated"] = display_df["last_updated"].apply(lambda x: x.strftime("%Y-%m-%d %H:%M:%S") if not pd.isna(x) else "")
    # Render with colored stage badges using markdown in dataframe is not directly supported; we'll render custom rows
    # Use an expander with a table and a per-row 'Inspect' action

    max_rows = st.number_input("Max rows to show (table)", min_value=5, max_value=500, value=200, step=5)
    show_df = display_df.head(max_rows)

    # Render as an interactive table using st.table but we want colored badges so we'll create an HTML table
    def render_html_table(df_to_render):
        header_cols = df_to_render.columns.tolist()
        html = ['<table style="width:100%; border-collapse:collapse;">']
        # header
        html.append("<thead><tr>")
        for c in header_cols:
            html.append(f'<th style="text-align:left;padding:6px;border-bottom:1px solid #ddd">{c}</th>')
        html.append("</tr></thead>")
        # rows
        html.append("<tbody>")
        for _, row in df_to_render.iterrows():
            html.append("<tr>")
            for c in header_cols:
                val = row[c]
                if c == "stage":
                    val_html = stage_color_html(val)
                elif c == "failed":
                    val_html = small_badge("failed", "#d90429") if val else small_badge("ok", "#1f9d55")
                elif c == "num_units_filtered":
                    if pd.isna(val) or val is None:
                        val_html = ""
                    else:
                        val_html = small_badge(str(val), "#f0ad4e")
                else:
                    val_html = f"<div style='font-size:0.95em'>{safe_str(val)}</div>"
                html.append(f'<td style="padding:6px;border-bottom:1px solid #f0f0f0;vertical-align:top">{val_html}</td>')
            html.append("</tr>")
        html.append("</tbody></table>")
        return "".join(html)

    st.markdown(render_html_table(show_df), unsafe_allow_html=True)

    st.markdown("---")
    st.subheader("Inspect a single checkpoint")

    # selection widget for single-checkpoint inspection
    all_files = view_df["file"].tolist()
    if not all_files:
        st.info("No checkpoints to inspect (after filters).")
        st.stop()

    sel_file = st.selectbox("Select checkpoint file", all_files, index=0)
    selected_row = view_df[view_df["file"] == sel_file].iloc[0]

    colA, colB = st.columns([3,2])
    with colA:
        st.markdown("**Metadata**")
        meta_display = {
            "file": selected_row["file"],
            "project": selected_row["project"],
            "date": selected_row["date"],
            "chip": selected_row["chip"],
            "run": selected_row["run"],
            "well": selected_row["well"],
            "rec": selected_row["rec"],
            "stage": selected_row["stage"],
            "failed": selected_row["failed"],
            "num_units_filtered": selected_row["num_units_filtered"],
            "last_updated": selected_row["last_updated"].strftime("%Y-%m-%d %H:%M:%S") if not pd.isna(selected_row["last_updated"]) else ""
        }
        st.json(meta_display)

        st.markdown("**Analyzer folder (click to copy)**")
        analyzer_folder = selected_row["_raw"].get("analyzer_folder") if isinstance(selected_row["_raw"], dict) else None
        if analyzer_folder:
            st.text_input("Analyzer folder", value=analyzer_folder, key="an_path")
            if st.button("Show in file browser (open externally)", key="open_fs"):
                st.info("Streamlit cannot open external file browsers. Copy the path above and paste into your file manager.")

        else:
            st.caption("No analyzer_folder field in JSON.")

    with colB:
        st.markdown("**Extensions computed**")
        exts = selected_row["_raw"].get("extensions_computed") if isinstance(selected_row["_raw"], dict) else []
        if exts:
            st.markdown(f"Computed extensions: **{len(exts)}**")
            # show simple bullet list
            for e in exts:
                st.write(f"- {e}")
        else:
            st.write("No extensions recorded.")

        st.markdown("**Error / Failure info**")
        if selected_row["failed"]:
            st.error(safe_str(selected_row["error"]))
            st.write(f"Failed flag: {selected_row['failed']}")
        else:
            st.success("No error recorded")

    st.markdown("**Full raw JSON**")
    st.expander("Raw JSON (expand)")  # small separator
    st.json(selected_row["_raw"])

    # Show parse/load errors if any
    if errors:
        st.markdown("---")
        st.error("Some checkpoint files failed to parse (see details below)")
        for fn, err in errors:
            st.write(f"- {fn}: {err}")

    st.markdown("---")
    st.caption("Tip: Use sidebar filters to narrow runs/wells. Use Export CSV to snapshot the current filtered view.")

# -------------------------
# Entrypoint
# -------------------------
if __name__ == "__main__":
    args = parse_args()
    cp_dir = args.checkpoint_dir if args.checkpoint_dir else str(Path(__file__).resolve().parent / "../AnalyzedData/checkpoints")
    run_app(cp_dir)