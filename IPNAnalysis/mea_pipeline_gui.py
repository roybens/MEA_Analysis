#!/usr/bin/env python3
"""
MEA Analysis Pipeline GUI
Launches run_pipeline_driver.py with selected arguments.
Place this file in the same directory as run_pipeline_driver.py
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import subprocess
import shlex
import os
from pathlib import Path

BASE_FILE_PATH = str(Path(__file__).resolve().parent)

# ── colour palette ────────────────────────────────────────────────────────────
BG        = "#0f1117"
PANEL     = "#1a1d27"
ACCENT    = "#4f8ef7"
ACCENT2   = "#2ecc71"
TEXT      = "#e8eaf0"
SUBTEXT   = "#7b8099"
BORDER    = "#2a2d3a"
DANGER    = "#e05252"
CHECK_ON  = "#4f8ef7"

class MEAGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MEA Analysis Pipeline")
        self.configure(bg=BG)
        self.resizable(True, True)

        # ── state vars ────────────────────────────────────────────────────────
        self.data_path      = tk.StringVar()
        self.output_dir     = tk.StringVar()
        self.checkpoint_dir = tk.StringVar()
        self.reference      = tk.StringVar()
        self.params         = tk.StringVar()
        self.docker         = tk.StringVar()
        self.sorter         = tk.StringVar(value="kilosort4")

        self.skip_sorting   = tk.BooleanVar()
        self.force_restart  = tk.BooleanVar()
        self.debug          = tk.BooleanVar()
        self.dry_run        = tk.BooleanVar()
        self.clean_up       = tk.BooleanVar()
        self.export_phy     = tk.BooleanVar()
        self.no_curation    = tk.BooleanVar()
        self.reanalyze      = tk.BooleanVar()
        self.fixed_y        = tk.BooleanVar()

        self._build_ui()

    # ── UI construction ───────────────────────────────────────────────────────
    def _build_ui(self):
        pad = {"padx": 20, "pady": 0}

        # Header
        hdr = tk.Frame(self, bg=ACCENT, height=4)
        hdr.pack(fill="x")

        title_frame = tk.Frame(self, bg=BG)
        title_frame.pack(fill="x", padx=20, pady=(18, 4))
        tk.Label(title_frame, text="MEA Analysis Pipeline",
                 font=("Courier New", 17, "bold"), bg=BG, fg=TEXT).pack(anchor="w")
        tk.Label(title_frame, text="Maxwell Biosystems · Spike Sorting & Burst Analysis",
                 font=("Courier New", 9), bg=BG, fg=SUBTEXT).pack(anchor="w")

        sep = tk.Frame(self, bg=BORDER, height=1)
        sep.pack(fill="x", padx=20, pady=(10, 0))

        # ── Path section ──────────────────────────────────────────────────────
        self._section("DATA PATH", pad)
        self._path_row("Data path (file or directory)*", self.data_path,
                       "dir", pad, required=True)
        self._path_row("Output directory*", self.output_dir,
                       "dir", pad, required=True)
        self._path_row("Checkpoint directory", self.checkpoint_dir, "dir", pad)
        self._path_row("Reference Excel file", self.reference, "file", pad)
        self._path_row("Params (JSON file or string)", self.params, "file", pad)

        sep2 = tk.Frame(self, bg=BORDER, height=1)
        sep2.pack(fill="x", padx=20, pady=(12, 0))

        # ── Text fields ───────────────────────────────────────────────────────
        self._section("SETTINGS", pad)
        text_frame = tk.Frame(self, bg=BG)
        text_frame.pack(fill="x", padx=20, pady=(0, 4))

        self._inline_entry(text_frame, "Sorter", self.sorter, col=0)
        self._inline_entry(text_frame, "Docker image", self.docker, col=2)

        sep3 = tk.Frame(self, bg=BORDER, height=1)
        sep3.pack(fill="x", padx=20, pady=(10, 0))

        # ── Checkboxes ────────────────────────────────────────────────────────
        self._section("OPTIONS", pad)

        checks = [
            (self.skip_sorting,  "Skip spike sorting",         "Run burst analysis only, skip Kilosort"),
            (self.reanalyze,     "Re-analyze bursts",          "Re-run burst analysis on existing spike times"),
            (self.fixed_y,       "Fixed Y axis",               "Use project y-max JSON to normalise all network plots"),
            (self.no_curation,   "No curation",                "Skip automatic unit curation"),
            (self.force_restart, "Force restart",              "Restart even if checkpoint says complete"),
            (self.export_phy,    "Export to Phy",              "Export sorted units to Phy format"),
            (self.clean_up,      "Clean up temp files",        "Delete binary / sorter output after run"),
            (self.dry_run,       "Dry run",                    "Print commands only, no processing"),
            (self.debug,         "Debug mode",                 "Verbose logging"),
        ]

        cb_frame = tk.Frame(self, bg=BG)
        cb_frame.pack(fill="x", padx=20, pady=(0, 6))

        for i, (var, label, tip) in enumerate(checks):
            row, col = divmod(i, 2)
            self._checkbox(cb_frame, var, label, tip, row=row, col=col * 2)

        sep4 = tk.Frame(self, bg=BORDER, height=1)
        sep4.pack(fill="x", padx=20, pady=(6, 0))

        # ── Command preview ───────────────────────────────────────────────────
        self._section("COMMAND PREVIEW", pad)
        self.cmd_text = tk.Text(self, height=3, bg=PANEL, fg=SUBTEXT,
                                font=("Courier New", 8), relief="flat",
                                wrap="word", state="disabled",
                                insertbackground=TEXT)
        self.cmd_text.pack(fill="x", padx=20, pady=(0, 8))

        # Update preview on any change
        for var in [self.data_path, self.output_dir, self.checkpoint_dir,
                    self.reference, self.params, self.docker, self.sorter,
                    self.skip_sorting, self.force_restart, self.debug,
                    self.dry_run, self.clean_up, self.export_phy,
                    self.no_curation, self.reanalyze, self.fixed_y]:
            var.trace_add("write", lambda *_: self._update_preview())

        # ── Run button ────────────────────────────────────────────────────────
        btn_frame = tk.Frame(self, bg=BG)
        btn_frame.pack(fill="x", padx=20, pady=(0, 20))

        run_btn = tk.Button(
            btn_frame, text="▶  RUN PIPELINE",
            font=("Courier New", 11, "bold"),
            bg=ACCENT, fg="white", activebackground="#3a7be0",
            activeforeground="white", relief="flat", cursor="hand2",
            padx=20, pady=10,
            command=self._run
        )
        run_btn.pack(side="left")

        copy_btn = tk.Button(
            btn_frame, text="Copy command",
            font=("Courier New", 9),
            bg=PANEL, fg=SUBTEXT, activebackground=BORDER,
            activeforeground=TEXT, relief="flat", cursor="hand2",
            padx=12, pady=10,
            command=self._copy_cmd
        )
        copy_btn.pack(side="left", padx=(8, 0))

        self._update_preview()

    # ── Helpers ───────────────────────────────────────────────────────────────
    def _section(self, title, pad):
        tk.Label(self, text=title, font=("Courier New", 8, "bold"),
                 bg=BG, fg=ACCENT).pack(anchor="w", padx=20, pady=(12, 4))

    def _path_row(self, label, var, kind, pad, required=False):
        row = tk.Frame(self, bg=BG)
        row.pack(fill="x", padx=20, pady=2)

        lbl_text = label + (" *" if required else "")
        tk.Label(row, text=lbl_text, width=30, anchor="w",
                 font=("Courier New", 9), bg=BG,
                 fg=TEXT if required else SUBTEXT).pack(side="left")

        entry = tk.Entry(row, textvariable=var, bg=PANEL, fg=TEXT,
                         insertbackground=TEXT, relief="flat",
                         font=("Courier New", 9), width=46)
        entry.pack(side="left", padx=(0, 6))

        def browse():
            if kind == "dir":
                p = filedialog.askdirectory()
            else:
                p = filedialog.askopenfilename()
            if p:
                var.set(p)

        tk.Button(row, text="Browse", font=("Courier New", 8),
                  bg=BORDER, fg=TEXT, activebackground=ACCENT,
                  activeforeground="white", relief="flat",
                  cursor="hand2", padx=8, pady=2,
                  command=browse).pack(side="left")

    def _inline_entry(self, parent, label, var, col):
        tk.Label(parent, text=label, font=("Courier New", 9),
                 bg=BG, fg=SUBTEXT, width=14, anchor="w").grid(
            row=0, column=col, sticky="w", padx=(0, 4))
        tk.Entry(parent, textvariable=var, bg=PANEL, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Courier New", 9), width=20).grid(
            row=0, column=col + 1, sticky="w", padx=(0, 20))

    def _checkbox(self, parent, var, label, tip, row, col):
        frm = tk.Frame(parent, bg=BG)
        frm.grid(row=row, column=col, sticky="w", padx=(0, 16), pady=3,
                 columnspan=2)

        cb = tk.Checkbutton(frm, variable=var, bg=BG,
                            activebackground=BG,
                            selectcolor=PANEL,
                            fg=CHECK_ON, activeforeground=CHECK_ON,
                            relief="flat", cursor="hand2")
        cb.pack(side="left")

        tk.Label(frm, text=label, font=("Courier New", 9, "bold"),
                 bg=BG, fg=TEXT).pack(side="left")
        tk.Label(frm, text=f"  —  {tip}", font=("Courier New", 8),
                 bg=BG, fg=SUBTEXT).pack(side="left")

    # ── Command builder ───────────────────────────────────────────────────────
    def _build_command(self):
        driver = os.path.join(BASE_FILE_PATH, "run_pipeline_driver.py")
        parts = ["python3", driver]

        path = self.data_path.get().strip()
        if path:
            parts.append(shlex.quote(path))

        if self.output_dir.get().strip():
            parts += ["--output-dir", shlex.quote(self.output_dir.get().strip())]
        if self.checkpoint_dir.get().strip():
            parts += ["--checkpoint-dir", shlex.quote(self.checkpoint_dir.get().strip())]
        if self.reference.get().strip():
            parts += ["--reference", shlex.quote(self.reference.get().strip())]
        if self.params.get().strip():
            parts += ["--params", shlex.quote(self.params.get().strip())]
        if self.docker.get().strip():
            parts += ["--docker", shlex.quote(self.docker.get().strip())]
        if self.sorter.get().strip() and self.sorter.get().strip() != "kilosort4":
            parts += ["--sorter", self.sorter.get().strip()]

        if self.skip_sorting.get():  parts.append("--skip-spikesorting")
        if self.force_restart.get(): parts.append("--force-restart")
        if self.debug.get():         parts.append("--debug")
        if self.dry_run.get():       parts.append("--dry")
        if self.clean_up.get():      parts.append("--clean-up")
        if self.export_phy.get():    parts.append("--export-to-phy")
        if self.no_curation.get():   parts.append("--no-curation")
        if self.reanalyze.get():     parts.append("--reanalyze-bursts")
        if self.fixed_y.get():       parts.append("--fixed-y")

        return " ".join(parts)

    def _update_preview(self):
        cmd = self._build_command()
        self.cmd_text.config(state="normal")
        self.cmd_text.delete("1.0", "end")
        self.cmd_text.insert("1.0", cmd)
        self.cmd_text.config(state="disabled")

    def _copy_cmd(self):
        self.clipboard_clear()
        self.clipboard_append(self._build_command())
        messagebox.showinfo("Copied", "Command copied to clipboard.")

    # ── Run ───────────────────────────────────────────────────────────────────
    def _run(self):
        if not self.data_path.get().strip():
            messagebox.showerror("Missing input", "Data path is required.")
            return
        if not self.output_dir.get().strip():
            messagebox.showerror("Missing input", "Output directory is required.")
            return

        cmd = self._build_command()
        print("\n" + "="*60)
        print("Launching MEA Pipeline:")
        print(cmd)
        print("="*60 + "\n")

        # Close GUI and hand control back to terminal
        self.destroy()
        os.system(cmd)


if __name__ == "__main__":
    app = MEAGui()
    app.mainloop()