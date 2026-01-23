#!/usr/bin/env python3

plates = set()

with open("raw_files.txt") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        # split S3 path
        parts = line.replace("s3://", "").split("/")

        # Expected structure:
        # 0 bucket
        # 1 raw
        # 2 experiment
        # 3 date
        # 4 plate
        if len(parts) < 5:
            continue

        bucket = parts[0]
        experiment = parts[2]
        date = parts[3]
        plate = parts[4]

        plate_path = f"s3://{bucket}/raw/{experiment}/{date}/{plate}/"
        plates.add(plate_path)

# Print sorted, one per line
for p in sorted(plates):
    print(p)