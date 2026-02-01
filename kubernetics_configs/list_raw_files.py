#!/usr/bin/env python3

import boto3

S3_ENDPOINT = "https://s3-west.nrp-nautilus.io"
BUCKET = "benshalom-lab-kcnt1-v5-primary"
PREFIX = ""   # optional: narrow path

s3 = boto3.client(
    "s3",
    endpoint_url=S3_ENDPOINT,
)

paginator = s3.get_paginator("list_objects_v2")

for page in paginator.paginate(Bucket=BUCKET, Prefix=PREFIX):
    for obj in page.get("Contents", []):
        key = obj["Key"]
        if key.endswith("data.raw.h5"):
            print(f"s3://{BUCKET}/{key}")