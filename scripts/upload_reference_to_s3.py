#!/usr/bin/env python3
"""
Upload reference data to S3 or Cloudflare R2 under methscope-data/reference/.
Run after building gene_list.json and gene_regions.db locally:
  python scripts/build_gene_regions.py
  python scripts/upload_reference_to_s3.py

For AWS S3: set S3_BUCKET and AWS credentials (AWS_ACCESS_KEY_ID, etc.).
For Cloudflare R2: set R2_BUCKET, R2_ACCOUNT_ID, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PREFIX = "methscope-data/reference"


def _get_client():
    import boto3
    from botocore.config import Config
    bucket = os.environ.get("R2_BUCKET") or os.environ.get("S3_BUCKET") or os.environ.get("AWS_S3_BUCKET")
    if not bucket:
        return None, None
    account_id = os.environ.get("R2_ACCOUNT_ID") or os.environ.get("CLOUDFLARE_ACCOUNT_ID")
    if account_id and (os.environ.get("R2_ACCESS_KEY_ID") or os.environ.get("R2_SECRET_ACCESS_KEY")):
        endpoint = f"https://{account_id}.r2.cloudflarestorage.com"
        client = boto3.client(
            "s3",
            endpoint_url=endpoint,
            region_name="auto",
            aws_access_key_id=os.environ.get("R2_ACCESS_KEY_ID", ""),
            aws_secret_access_key=os.environ.get("R2_SECRET_ACCESS_KEY", ""),
            config=Config(signature_version="s3v4"),
        )
        return client, bucket
    client = boto3.client("s3", region_name=os.environ.get("AWS_REGION", "us-east-1"))
    return client, bucket


def main():
    try:
        import boto3
    except ImportError:
        print("Install boto3: pip install boto3", file=sys.stderr)
        sys.exit(1)
    client, bucket = _get_client()
    if not client or not bucket:
        print("Set R2_BUCKET + R2_ACCOUNT_ID + R2_ACCESS_KEY_ID + R2_SECRET_ACCESS_KEY", file=sys.stderr)
        print("  or  S3_BUCKET + AWS credentials.", file=sys.stderr)
        sys.exit(1)

    to_upload = [
        ("gene_list.json", "gene_list.json"),
        ("gene_regions.db", "gene_regions.db"),
        ("genomic.gff", "genomic.gff"),
    ]
    for local_name, s3_name in to_upload:
        path = DATA / local_name
        if not path.exists():
            print(f"Skip (not found): {path}")
            continue
        key = f"{PREFIX}/{s3_name}"
        print(f"Uploading {path} -> s3://{bucket}/{key}")
        client.upload_file(str(path), bucket, key)
    print("Done.")


if __name__ == "__main__":
    main()
