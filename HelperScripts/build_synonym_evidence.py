#!/usr/bin/env python3
"""Build the FlyBase symbol-synonym evidence index from the Chado XML gene file.

For every Dmel gene, the index lists each symbol synonym together with the number of
distinct FlyBase publications (FBrf) that used it. GetFBgnIDs.py uses these counts to
resolve a synonym shared by several genes without calling any external service.

Usage:
    python HelperScripts/build_synonym_evidence.py              # download if needed, then build
    python HelperScripts/build_synonym_evidence.py --skip-download
    python HelperScripts/build_synonym_evidence.py --sample path/to/prefix.xml.gz --output out.tsv.gz

The XML is machine-generated with a fixed layout (see chado_FBgn.dtd), so it is read
as a stream of lines with a small state machine instead of building an XML tree.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable, Iterator

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

FLYBASE_RELEASE = "FB2026_03"
CHADO_XML_BASE_URL = f"https://s3ftp.flybase.org/releases/{FLYBASE_RELEASE}/chado-xml"
CHADO_DIR = PROJECT_ROOT / "Data" / "FlyBase_Chado"
CHADO_XML_PATH = CHADO_DIR / "chado_FBgn.xml.gz"
EVIDENCE_PATH = CHADO_DIR / f"synonym_evidence_{FLYBASE_RELEASE}.tsv.gz"
EVIDENCE_COLUMNS = ["fbgn", "current_symbol", "synonym", "n_pubs", "is_current_symbol"]

_XML_ENTITIES = (("&lt;", "<"), ("&gt;", ">"), ("&quot;", '"'), ("&apos;", "'"), ("&amp;", "&"))


def _unescape(text: str) -> str:
    for entity, char in _XML_ENTITIES:
        text = text.replace(entity, char)
    return text


def _tag_value(line: bytes, prefix: bytes, closing: bytes) -> str:
    return _unescape(line[len(prefix):line.rindex(closing)].decode("utf-8", "replace"))


def parse_genes(lines: Iterable[bytes]) -> Iterator[tuple[str, str, dict[str, int]]]:
    """Yield (fbgn, current_symbol, {symbol_synonym: n_distinct_pubs}) for each Dmel gene.

    Top-level features sit one tab deep; their direct fields two tabs deep. Inside a
    feature_synonym block, synonym fields are five tabs deep, and the synonym/pub
    wrappers are bracketed by four-tab <synonym>/<pub> lines.
    """
    in_feature = False
    uniquename = name = organism = ""
    obsolete = False
    in_syn_block = in_synonym = in_pub = False
    syn_is_symbol = False
    syn_name = pub_id = ""
    is_internal = False
    is_current = False
    pubs: dict[str, set[str]] = {}
    current_symbol = ""

    for line in lines:
        if not in_feature:
            if line.startswith(b"\t<feature>"):
                in_feature = True
                uniquename = name = organism = current_symbol = ""
                obsolete = False
                in_syn_block = in_synonym = in_pub = False
                pubs = {}
            continue
        if line.startswith(b"\t</feature>"):
            in_feature = False
            if uniquename.startswith("FBgn") and organism == "Dmel" and not obsolete:
                yield uniquename, current_symbol or name, {syn: len(ids) for syn, ids in pubs.items()}
            continue
        if not line.startswith(b"\t\t"):
            continue
        if line.startswith(b"\t\t<uniquename>"):
            uniquename = _tag_value(line, b"\t\t<uniquename>", b"</uniquename>")
        elif line.startswith(b"\t\t<name>"):
            name = _tag_value(line, b"\t\t<name>", b"</name>")
        elif line.startswith(b"\t\t<is_obsolete>"):
            obsolete = line.startswith(b"\t\t<is_obsolete>1")
        elif line.startswith(b"\t\t\t\t<abbreviation>"):
            organism = _tag_value(line, b"\t\t\t\t<abbreviation>", b"</abbreviation>")
        elif line.startswith(b"\t\t<feature_synonym>"):
            in_syn_block = True
            in_synonym = in_pub = syn_is_symbol = False
            syn_name = pub_id = ""
            is_internal = is_current = False
        elif in_syn_block:
            if line.startswith(b"\t\t</feature_synonym>"):
                in_syn_block = False
                if syn_is_symbol and syn_name:
                    ids = pubs.setdefault(syn_name, set())
                    if pub_id and not is_internal:
                        ids.add(pub_id)
                    if is_current:
                        current_symbol = syn_name
            elif line.startswith(b"\t\t\t\t<synonym>"):
                in_synonym = True
            elif line.startswith(b"\t\t\t\t</synonym>"):
                in_synonym = False
            elif line.startswith(b"\t\t\t\t<pub>"):
                in_pub = True
            elif line.startswith(b"\t\t\t\t</pub>"):
                in_pub = False
            elif in_synonym and line.startswith(b"\t\t\t\t\t\t\t<name>"):
                syn_is_symbol = line.startswith(b"\t\t\t\t\t\t\t<name>symbol</name>")
            elif in_synonym and line.startswith(b"\t\t\t\t\t<name>"):
                syn_name = _tag_value(line, b"\t\t\t\t\t<name>", b"</name>")
            elif in_pub and line.startswith(b"\t\t\t\t\t<uniquename>"):
                pub_id = _tag_value(line, b"\t\t\t\t\t<uniquename>", b"</uniquename>")
            elif line.startswith(b"\t\t\t<is_internal>"):
                is_internal = line.startswith(b"\t\t\t<is_internal>1")
            elif line.startswith(b"\t\t\t<is_current>"):
                is_current = line.startswith(b"\t\t\t<is_current>1")


def _gunzip_lines(path: Path) -> Iterator[bytes]:
    """Stream decompressed lines using the system gzip (faster than gzip.open)."""
    proc = subprocess.Popen(["gzip", "-dc", str(path)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, bufsize=1 << 22)
    try:
        assert proc.stdout is not None
        yield from proc.stdout
    finally:
        if proc.stdout:
            proc.stdout.close()
        proc.wait()


def build_index(xml_path: Path, output_path: Path) -> int:
    """Stream the gene XML and write the evidence index. Returns the number of genes."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.tmp")
    started = time.time()
    genes = rows = 0
    try:
        with gzip.open(tmp_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(EVIDENCE_COLUMNS)
            for fbgn, current, synonyms in parse_genes(_gunzip_lines(xml_path)):
                genes += 1
                for synonym, n_pubs in sorted(synonyms.items()):
                    writer.writerow([fbgn, current, synonym, n_pubs, int(synonym == current)])
                    rows += 1
                if genes % 2000 == 0:
                    print(f"  {genes:,} genes, {rows:,} synonym rows, {time.time() - started:,.0f}s", flush=True)
        os.replace(tmp_path, output_path)
    finally:
        tmp_path.unlink(missing_ok=True)
    print(f"Wrote {output_path}: {genes:,} genes, {rows:,} synonym rows in {time.time() - started:,.0f}s")
    return genes


def download(url: str, destination: Path):
    """Resumable download (HTTP Range) of a large file."""
    import requests

    destination.parent.mkdir(parents=True, exist_ok=True)
    have = destination.stat().st_size if destination.exists() else 0
    headers = {"Range": f"bytes={have}-"} if have else {}
    with requests.get(url, headers=headers, stream=True, timeout=60) as response:
        if response.status_code == 416:
            return
        response.raise_for_status()
        mode = "ab" if response.status_code == 206 else "wb"
        with open(destination, mode) as handle:
            for chunk in response.iter_content(1 << 20):
                handle.write(chunk)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sample", help="Parse this local (possibly truncated) .xml.gz instead of the full file.")
    parser.add_argument("--output", help="Index path (default: Data/FlyBase_Chado/synonym_evidence_<release>.tsv.gz).")
    parser.add_argument("--skip-download", action="store_true", help="Do not download the Chado XML if it is missing.")
    args = parser.parse_args(argv)

    xml_path = Path(args.sample) if args.sample else CHADO_XML_PATH
    output = Path(args.output) if args.output else EVIDENCE_PATH
    if not args.sample and not args.skip_download:
        for name in ("chado_FBgn.dtd", "chado_FBgn.xml.gz"):
            target = CHADO_DIR / name
            print(f"Ensuring {target.name} (resumable download)")
            download(f"{CHADO_XML_BASE_URL}/{name}", target)
    if not xml_path.exists():
        print(f"Error: Chado XML not found: {xml_path}")
        return 1
    if not args.sample:
        print("Checking gzip integrity ...")
        if subprocess.run(["gzip", "-t", str(xml_path)]).returncode != 0:
            print("Error: chado_FBgn.xml.gz is corrupt or incomplete; rerun to resume the download.")
            return 1
    build_index(xml_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
