from __future__ import annotations

import argparse
from pathlib import Path

from .reporting import UnitMatchReportConfig, generate_unitmatch_static_report_pack


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate static UnitMatch report pack from existing artifacts",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=str,
        help="Well-level output directory containing UnitMatch artifacts",
    )
    parser.add_argument(
        "--throughput-subdir-name",
        default="unitmatch_throughput",
        type=str,
        help="Throughput subdirectory name (default: unitmatch_throughput)",
    )
    parser.add_argument(
        "--output-subdir-name",
        default="unitmatch_outputs",
        type=str,
        help="UnitMatch output subdirectory name (default: unitmatch_outputs)",
    )
    parser.add_argument(
        "--report-subdir-name",
        default="unitmatch_reports",
        type=str,
        help="Report output subdirectory name (default: unitmatch_reports)",
    )
    parser.add_argument(
        "--max-heatmap-units",
        default=200,
        type=int,
        help="Maximum units displayed in similarity heatmap (default: 200)",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    cfg = UnitMatchReportConfig(
        throughput_subdir_name=str(args.throughput_subdir_name),
        output_subdir_name=str(args.output_subdir_name),
        report_subdir_name=str(args.report_subdir_name),
        max_heatmap_units=int(args.max_heatmap_units),
    )

    report = generate_unitmatch_static_report_pack(
        output_dir=output_dir,
        logger=None,
        config=cfg,
    )
    print("UnitMatch report generation complete")
    print(f"Report directory: {report.get('report_root')}")


if __name__ == "__main__":
    main()
