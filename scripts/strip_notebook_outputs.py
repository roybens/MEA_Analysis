#!/usr/bin/env python3
"""Remove cell outputs and execution counts from Jupyter notebooks."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable


def iter_notebooks(argv: list[str]) -> Iterable[Path]:
    if argv:
        for name in argv:
            path = Path(name)
            if path.suffix == ".ipynb" and path.exists():
                yield path
        return

    yield from Path(".").rglob("*.ipynb")


def strip_notebook(path: Path) -> bool:
    original_text = path.read_text(encoding="utf-8")
    notebook = json.loads(original_text)

    changed = False
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue

        if cell.get("outputs"):
            cell["outputs"] = []
            changed = True

        if cell.get("execution_count") is not None:
            cell["execution_count"] = None
            changed = True

    metadata = notebook.get("metadata", {})
    if "widgets" in metadata:
        metadata.pop("widgets", None)
        changed = True

    if not changed:
        return False

    path.write_text(
        json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8",
    )
    return True


def main() -> int:
    changed_files: list[str] = []

    for notebook_path in iter_notebooks(sys.argv[1:]):
        if strip_notebook(notebook_path):
            changed_files.append(str(notebook_path))

    if changed_files:
        print("Stripped notebook outputs:")
        for changed_file in changed_files:
            print(f"  {changed_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
