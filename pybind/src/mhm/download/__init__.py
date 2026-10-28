r"""Download runnable v6 example domains from the mHM repository.

\ingroup py_mhm
"""

import argparse
import shutil
import tarfile
from gzip import BadGzipFile
from pathlib import Path, PurePosixPath, PureWindowsPath
from tempfile import TemporaryDirectory
from urllib.error import URLError
from urllib.parse import quote, urlencode
from urllib.request import urlretrieve

MHM_URLS = {
    "github": "https://github.com/mhm-ufz/mHM/archive/{branch}.tar.gz",
    "git.ufz": "https://git.ufz.de/mhm/mhm/-/archive/{branch}/mhm-{branch}.tar.gz?{query}",
}
VALID_DOMAINS = [1, 2]
VALID_SOURCES = ["auto", "github", "git.ufz"]
SOURCE_ORDER = {
    "auto": ("github", "git.ufz"),
    "github": ("github",),
    "git.ufz": ("git.ufz",),
}


class _InvalidArchiveError(ValueError):
    """A repository archive does not contain a usable example domain."""


def _extract_domain(archive, extracted, name):
    """Extract only regular files/directories beneath a unique example root."""
    found = False
    roots = set()
    extracted.mkdir()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            member_path = PurePosixPath(member.name)
            if (
                member_path.is_absolute()
                or ".." in member_path.parts
                or "\\" in member.name
            ):
                raise _InvalidArchiveError(
                    f"Invalid path in example archive: {member.name}"
                )
            parts = member_path.parts
            if len(parts) < 3 or parts[1:3] != ("examples", name):
                continue
            if len(parts) == 3 and not member.isdir():
                raise _InvalidArchiveError("Example root is not a directory")
            if PureWindowsPath(*parts[3:]).drive:
                raise _InvalidArchiveError(
                    f"Invalid path in example archive: {member.name}"
                )
            roots.add(parts[0])
            if len(roots) > 1:
                raise _InvalidArchiveError("Archive contains multiple example roots")
            destination = extracted.joinpath(*parts[3:])
            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                destination.parent.mkdir(parents=True, exist_ok=True)
                with tar.extractfile(member) as source:
                    with destination.open("wb") as target:
                        shutil.copyfileobj(source, target)
                found = True
            else:
                raise _InvalidArchiveError(
                    f"Unsupported entry in example archive: {member.name}"
                )
    if not found:
        raise _InvalidArchiveError(f"Archive does not contain examples/{name}")


def download_examples(branch=None, domain=1, path=None, verbose=False, source="auto"):
    """Download one example domain, preserving its runnable directory layout.

    Development packages use v6alpha; releases use their matching version tag.
    Explicit branch/tag/commit and destination arguments override these defaults.
    Source "auto" tries GitHub first, then git.ufz; either source can be selected
    explicitly. Local destination errors do not trigger another download.
    """
    from .. import __version__

    if domain not in VALID_DOMAINS:
        raise ValueError(f"mhm-download: 'domain' needs to be 1 or 2. Got: '{domain}'")
    if source not in VALID_SOURCES:
        raise ValueError(
            f"mhm-download: 'source' needs to be one of {VALID_SOURCES}. Got: '{source}'"
        )
    if branch is None:
        branch = "v6alpha" if "dev" in __version__ else f"v{__version__}"
    name = f"domain_{domain:02d}"
    folder = f"examples/{name}"
    path = Path(path or name)
    if verbose:
        print(
            f"Downloading mHM example domain {domain} from {branch} to {path.resolve()}"
        )
    failures = []
    with TemporaryDirectory() as tmp_dir:
        for index, current_source in enumerate(SOURCE_ORDER[source]):
            archive = Path(tmp_dir) / f"example-{index}.tar.gz"
            extracted = Path(tmp_dir) / f"example-{index}"
            url = MHM_URLS[current_source].format(
                branch=quote(branch, safe=""), query=urlencode({"path": folder})
            )
            if verbose:
                print(f"  trying source: {current_source}")
            try:
                urlretrieve(url, archive)
                _extract_domain(archive, extracted, name)
            except (
                URLError,
                tarfile.TarError,
                BadGzipFile,
                EOFError,
                _InvalidArchiveError,
            ) as error:
                failures.append((current_source, error))
                if verbose:
                    print(f"  source failed: {current_source} ({error})")
                continue
            # Copy errors concern the destination, not the remote archive.
            if verbose:
                print(f"  selected source: {current_source}")
            shutil.copytree(extracted, path)
            return

    attempts = "; ".join(
        f"{name}: {type(error).__name__}: {error}" for name, error in failures
    )
    raise RuntimeError(
        f"mhm-download: failed to download {folder} at ref '{branch}'. Attempts: {attempts}"
    ) from failures[-1][1]


def cli(argv=None):
    """Command line interface for example-domain downloads."""
    from .. import __version__

    parser = argparse.ArgumentParser(
        prog="mhm-download",
        description="Download a runnable mHM v6 example domain.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-V", "--version", action="version", version=__version__)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-b", "--branch", help="Repository branch, tag, or commit")
    parser.add_argument("-d", "--domain", type=int, default=1, choices=VALID_DOMAINS)
    parser.add_argument(
        "-p", "--path", help="Destination; defaults to domain_01 or domain_02"
    )
    parser.add_argument(
        "-s",
        "--source",
        choices=VALID_SOURCES,
        default="auto",
        help="Download source; auto tries GitHub first, then git.ufz",
    )
    args = parser.parse_args(argv)
    download_examples(args.branch, args.domain, args.path, args.verbose, args.source)
