from pathlib import Path
import hashlib

template = Path("mea_spikesorter_job_template.yaml").read_text()
paths = Path("plate_paths.txt").read_text().splitlines()

outdir = Path("jobs")
outdir.mkdir(exist_ok=True)

for s3 in paths:
    s3 = s3.strip().rstrip("/")
    if not s3:
        continue

    parts = s3.split("/")
    date = parts[-2]
    plate = parts[-1]

    job_id = f"{plate.lower()}-{date}"

    yaml = (
        template
        .replace("{{ JOB_ID }}", job_id)
        .replace("{{ PLATE_ID }}", plate)
        .replace("{{ DATE }}", date)
        .replace("{{ S3_INPUT_PATH }}", s3 + "/")
        .replace("{{ S3_OUTPUT_PATH}}", s3.replace("/raw/", "/processed/") + "/")
    )

    outfile = outdir / f"job_{job_id}.yaml"
    outfile.write_text(yaml)

print(f"Generated {len(list(outdir.glob('*.yaml')))} jobs")