"""FlyBase data path resolution and bootstrap utilities."""

from __future__ import annotations

import os
from glob import glob
from pathlib import Path

import requests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOCAL_DATA_ROOT = PROJECT_ROOT / "Data"

_TURBO_CANDIDATES = [
    Path("/Volumes/umms-rallada"),
    Path("/nfs/turbo/umms-rallada"),
]
TURBO_ROOT = next((p for p in _TURBO_CANDIDATES if p.exists()), _TURBO_CANDIDATES[0])
SHARED_DATA_ROOT = TURBO_ROOT / "UM Lab Users" / "Aadish" / "Data"

FLYBASE_REQUIRED_FILES = [
    {
        "label": "FlyBase synonyms",
        "subdir": "Genes",
        "pattern": "fb_synonym",
        "url": "https://s3ftp.flybase.org/releases/current/precomputed_files/synonyms/fb_synonym_fb_2026_01.tsv.gz",
    },
    {
        "label": "FlyBase entity-publication associations",
        "subdir": "FlyBase_References",
        "pattern": "entity_publication",
        "url": "https://s3ftp.flybase.org/releases/current/precomputed_files/references/entity_publication_fb_2026_01.tsv.gz",
    },
    {
        "label": "FlyBase PMID/PMCID/DOI mapping",
        "subdir": "FlyBase_References",
        "pattern": "fbrf_pmid_pmcid_doi",
        "url": "https://s3ftp.flybase.org/releases/current/precomputed_files/references/fbrf_pmid_pmcid_doi_fb_2026_01.tsv.gz",
    },
]


def _is_accessible_dir(path: Path) -> bool:
    try:
        return path.exists() and path.is_dir() and os.access(path, os.R_OK | os.X_OK)
    except OSError:
        return False


def default_data_root() -> Path:
    """Use shared lab data when available, otherwise repo-local Data."""
    env_root = os.environ.get("FLAI_DATA_DIR")
    if env_root:
        return Path(env_root).expanduser()
    if _is_accessible_dir(SHARED_DATA_ROOT):
        return SHARED_DATA_ROOT
    return LOCAL_DATA_ROOT


def using_local_data_root() -> bool:
    return default_data_root().resolve() == LOCAL_DATA_ROOT.resolve()


def resolve_flybase_data_dir(explicit: str | os.PathLike[str] | None = None) -> Path:
    if explicit:
        return Path(explicit).expanduser()
    env_path = os.environ.get("FLYBASE_DATA_DIR")
    if env_path:
        return Path(env_path).expanduser()
    return default_data_root() / "FlyBase"


def resolve_pubmed_cache_dir(explicit: str | os.PathLike[str] | None = None) -> Path:
    if explicit:
        return Path(explicit).expanduser()
    env_path = os.environ.get("PUBMED_CACHE_DIR") or os.environ.get("FLAI_PUBMED_CACHE_DIR")
    if env_path:
        return Path(env_path).expanduser()
    return default_data_root() / "PubMed Cache"


def _has_matching_tsv(directory: Path, pattern: str) -> bool:
    return bool(glob(str(directory / f"{pattern}*.tsv.gz")) or glob(str(directory / f"{pattern}*.tsv")))


def missing_flybase_files(flybase_data_dir: Path) -> list[dict[str, str]]:
    missing = []
    for requirement in FLYBASE_REQUIRED_FILES:
        directory = flybase_data_dir / requirement["subdir"]
        if not _has_matching_tsv(directory, requirement["pattern"]):
            missing.append(requirement)
    return missing


def _download_file(url: str, dest: Path, timeout: int = 300) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with open(tmp, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp.replace(dest)


def ensure_flybase_data_files(flybase_data_dir: Path, download_missing: bool = True) -> list[dict[str, str]]:
    """Download missing FlyBase TSVs when possible and return anything still missing."""
    flybase_data_dir = Path(flybase_data_dir)
    missing = missing_flybase_files(flybase_data_dir)
    if not missing or not download_missing:
        return missing

    print(f"FlyBase data missing under: {flybase_data_dir}")
    print("Attempting to download current FlyBase TSVs into the local expected layout...")

    for requirement in missing:
        url = requirement["url"]
        dest = flybase_data_dir / requirement["subdir"] / Path(url).name
        try:
            print(f"  Downloading {requirement['label']}: {dest.name}")
            _download_file(url, dest)
        except Exception as exc:
            print(f"  [Warning] Could not download {requirement['label']} from {url}: {exc}")

    return missing_flybase_files(flybase_data_dir)
