#!/usr/bin/env python3
"""Build the single-file Drug Readout page from editable HTML and local CSV data.

Like the Reverse Metabolomics browser build, the editable page has empty data
blocks and the distributed page embeds compressed data in those blocks. Run:

    python3 webapp/build_drug_readouts.py
    python3 webapp/build_drug_readouts.py --check
    python3 webapp/build_drug_readouts.py --zenodo

The first build of this port used --bootstrap once to extract an editable
source from the existing generated page. Normal updates should edit the source
HTML or data/*.csv, then rerun the builder. Only Python's standard library is
needed. gzip mtime=0 makes unchanged inputs produce byte-identical output.
"""

import argparse
import base64
import csv
import gzip
import hashlib
import io
import os
import re
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "source" / "drug_readouts.html"
OUTPUT = HERE / "drug_readouts.html"
DATA = (
    ("drug-metadata", ROOT / "data" / "GNPS_Drug_Library_Metadata_Drugs.csv",
     {"gnps_libid", "name_parent_compound", "chemical_source",
      "pharmacologic_class", "therapeutic_area", "therapeutic_indication"}),
    ("analog-metadata", ROOT / "data" / "GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv",
     {"analog_libid", "name_connected_compound", "chemical_source",
      "pharmacologic_class", "therapeutic_area", "therapeutic_indication"}),
)
ZENODO_RECORD = "https://zenodo.org/records/17232042"
ZENODO_FILES = (
    ("GNPS_Drug_Library_Metadata_Drugs.csv", "1770aec351079484e66f604c1fb1176f"),
    ("GNPS_Drug_Library_Metadata_Drug_Analogs_Updated_With_Source.csv",
     "16a31645a5c2334bbdc8a9734c63e62f"),
)


def fail(message):
    raise SystemExit("build_drug_readouts: " + message)


def block_pattern(block_id):
    return re.compile(
        r'(<script type="application/octet-stream" id="' + re.escape(block_id)
        + r'"[^>]*>)(.*?)(</script>)', re.DOTALL
    )


def replace_block(page, block_id, replacement, require_empty=False):
    pattern = block_pattern(block_id)
    matches = list(pattern.finditer(page))
    if len(matches) != 1:
        fail(f'page needs exactly one <script id="{block_id}"> block')
    if require_empty and matches[0].group(2).strip():
        fail(f'{block_id} is not empty in {SOURCE}; edit the source, not the generated page')
    return pattern.sub(lambda _match: replacement, page, count=1)


def bootstrap_source():
    if SOURCE.exists():
        fail(f"refusing to overwrite existing editable source: {SOURCE}")
    page = OUTPUT.read_text(encoding="utf-8")
    for block_id, _path, _columns in DATA:
        empty = f'<script type="application/octet-stream" id="{block_id}">\n</script>'
        page = replace_block(page, block_id, empty)
    SOURCE.parent.mkdir(parents=True, exist_ok=True)
    with SOURCE.open("x", encoding="utf-8") as stream:
        stream.write(page)
    print(f"build_drug_readouts: created editable source {SOURCE}")


def compressed_block(block_id, path, required_columns):
    if not path.is_file():
        fail(f"data file not found: {path}")
    raw = path.read_bytes()
    if not raw:
        fail(f"data file is empty: {path}")
    header = next(csv.reader([raw.splitlines()[0].decode("utf-8-sig")]))
    missing = required_columns.difference(header)
    if missing:
        fail(f"{path.name} is missing required column(s): {', '.join(sorted(missing))}")
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", filename="", compresslevel=9, mtime=0) as stream:
        stream.write(raw)
    encoded = base64.encodebytes(buffer.getvalue()).decode("ascii")
    digest = hashlib.sha256(raw).hexdigest()
    block = (
        f'<script type="application/octet-stream" id="{block_id}"'
        f' data-source="{path.name}" data-sha256="{digest}">\n'
        f'{encoded}</script>'
    )
    return block, len(raw), len(buffer.getvalue())


def download_zenodo(directory):
    paths = []
    for filename, expected_md5 in ZENODO_FILES:
        url = f"{ZENODO_RECORD}/files/{filename}?download=1"
        target = directory / filename
        print(f"build_drug_readouts: downloading {filename} from Zenodo…", file=sys.stderr)
        try:
            with urllib.request.urlopen(url, timeout=120) as response, target.open("wb") as stream:
                shutil.copyfileobj(response, stream)
        except (OSError, urllib.error.URLError) as error:
            fail(f"could not download {filename}: {error}")
        actual_md5 = hashlib.md5(target.read_bytes()).hexdigest()
        if actual_md5 != expected_md5:
            fail(f"Zenodo checksum mismatch for {filename}: expected {expected_md5}, got {actual_md5}")
        paths.append(target)
    return paths


def build(paths):
    page = SOURCE.read_text(encoding="utf-8")
    sizes = []
    for (block_id, _default, columns), path in zip(DATA, paths):
        block, raw_size, zipped_size = compressed_block(block_id, path, columns)
        page = replace_block(page, block_id, block, require_empty=True)
        sizes.append((path.name, raw_size, zipped_size))
    for token in ("{{", "{%", "{#"):
        if token in page:
            fail(f"generated HTML contains Jinja token {token!r}")
    return page.encode("utf-8"), sizes


def write_atomically(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_name = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=path.parent, prefix=".drug_readouts-",
                                         suffix=".tmp", delete=False) as stream:
            temp_name = stream.name
            stream.write(data)
        # NamedTemporaryFile creates 0600; web servers need to read the page.
        os.chmod(temp_name, 0o644)
        os.replace(temp_name, path)
    finally:
        if temp_name and os.path.exists(temp_name):
            os.unlink(temp_name)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="fail if the generated page is missing or stale")
    parser.add_argument("--bootstrap", action="store_true", help="one-time extraction of editable HTML from the existing page")
    parser.add_argument("--output", type=Path, default=OUTPUT, help="generated HTML destination")
    parser.add_argument("--drugs", type=Path, default=DATA[0][1], help="drug metadata CSV (default: repository data file)")
    parser.add_argument("--analogs", type=Path, default=DATA[1][1], help="drug-analog metadata CSV (default: repository data file)")
    parser.add_argument("--zenodo", action="store_true", help="download and embed the two metadata CSVs from Zenodo record 17232042")
    args = parser.parse_args()
    if args.bootstrap:
        if args.check or args.zenodo or args.output != OUTPUT or args.drugs != DATA[0][1] or args.analogs != DATA[1][1]:
            fail("--bootstrap cannot be combined with other options")
        bootstrap_source()
        return
    if args.zenodo and (args.drugs != DATA[0][1] or args.analogs != DATA[1][1]):
        fail("--zenodo cannot be combined with --drugs or --analogs")
    if args.zenodo:
        with tempfile.TemporaryDirectory(prefix="drug-readouts-zenodo-") as directory:
            page, sizes = build(download_zenodo(Path(directory)))
    else:
        page, sizes = build((args.drugs, args.analogs))
    if args.check:
        if not args.output.exists() or args.output.read_bytes() != page:
            fail(f"{args.output} is missing or stale; rerun this script without --check")
        print(f"build_drug_readouts: {args.output} is up to date")
        return
    if args.output.exists() and args.output.read_bytes() == page:
        print(f"build_drug_readouts: {args.output} is already up to date")
        return
    write_atomically(args.output, page)
    detail = ", ".join(f"{name}: {raw_size // 1024} KiB -> {zipped_size // 1024} KiB"
                       for name, raw_size, zipped_size in sizes)
    print(f"build_drug_readouts: wrote {args.output} ({len(page) // 1024} KiB; {detail})")


if __name__ == "__main__":
    main()
