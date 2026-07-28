"""!
Download routines to get test domains for mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
@ingroup mhm
"""

import argparse
import shutil
import tarfile
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from urllib.error import ContentTooShortError, HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlretrieve

MHM_URLS = {
    "git.ufz": (
        "https://git.ufz.de/mhm/mhm/-/archive/{branch}/"
        "mhm-{branch}.{format}?path={folder}"
    ),
    "github": "https://github.com/mhm-ufz/mHM/archive/{branch}.{format}",
}
VALID_DOMAINS = [1, 2]
VALID_SOURCES = ["auto", "git.ufz", "github"]
SOURCE_ORDER = {
    "auto": ("git.ufz", "github"),
    "git.ufz": ("git.ufz",),
    "github": ("github",),
}


class _InvalidArchiveError(Exception):
    """An archive does not contain the requested test domain."""


def _url(source, branch, format, folder):
    branch = quote(branch, safe="")
    return MHM_URLS[source].format(
        branch=branch,
        format=format,
        folder=folder,
    )


def _dl(source, branch, format, folder, filename):
    urlretrieve(_url(source, branch, format, folder), filename)


def _extract_domain(tar_file, tar_dir, folder):
    """Extract and return a test-domain directory from a repository archive."""
    with tarfile.open(tar_file, "r:gz") as tar:
        members = []
        roots = set()
        for member in tar.getmembers():
            member_path = PurePosixPath(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise _InvalidArchiveError(
                    f"archive contains unsafe path: '{member.name}'"
                )
            if len(member_path.parts) < 2 or member_path.parts[1] != folder:
                continue
            roots.add(member_path.parts[0])
            members.append(member)

        if not members or len(roots) != 1:
            raise _InvalidArchiveError(
                f"archive does not contain a unique '{folder}' directory"
            )

        filter_args = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}
        tar.extractall(path=tar_dir, members=members, **filter_args)

    folder_path = tar_dir / roots.pop() / folder
    if not folder_path.is_dir():
        raise _InvalidArchiveError(
            f"archive does not contain a usable '{folder}' directory"
        )
    return folder_path


def download_test(branch=None, domain=1, path=None, verbose=False, source="auto"):
    """
    Download a test domain for mHM.

    @param branch (str, optional): Branch, tag, or commit of the mHM repository to
        take the test domain from, by default tag determined from the mHM version
    @param domain (int, optional): Test domain 1 or 2, by default 1
    @param path (pathlike, optional): Destination path for the downloaded folder,
        by default original name of the test domain folder
    @param verbose (bool, optional): Report download details, by default False
    @param source (str, optional): Download source: "auto", "git.ufz", or "github";
        "auto" tries git.ufz first and then github, by default "auto"
    """
    # format fixed to tar.gz
    format = "tar.gz"
    # determine branch from mhm version
    if branch is None:
        from .. import __version__

        branch = "main" if "dev" in __version__ else f"v{__version__}"
    # check test domain and source
    if domain not in VALID_DOMAINS:
        msg = f"mhm-download: 'domain' needs to be 1 or 2. Got: '{domain}'"
        raise ValueError(msg)
    if source not in VALID_SOURCES:
        msg = (
            "mhm-download: 'source' needs to be 'auto', 'git.ufz', or 'github'. "
            f"Got: '{source}'"
        )
        raise ValueError(msg)

    folder = "test_domain" if domain == 1 else "test_domain_2"
    path = Path(path or folder)
    if verbose:
        print(f"Downloading mHM test domain '{domain}'")
        print(f"  branch: '{branch}'")
        print(f"  source: '{source}'")
        print(f"  path:   '{path.resolve()}'")

    failures = []
    source_errors = (
        ContentTooShortError,
        HTTPError,
        URLError,
        tarfile.TarError,
        _InvalidArchiveError,
    )
    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        for index, current_source in enumerate(SOURCE_ORDER[source]):
            tar_file = tmp_dir / f"test-{index}.{format}"
            tar_dir = tmp_dir / f"test-{index}"
            if verbose:
                print(f"  trying source: '{current_source}'")
            try:
                _dl(current_source, branch, format, folder, tar_file)
                folder_path = _extract_domain(tar_file, tar_dir, folder)
            except source_errors as err:
                failures.append((current_source, err))
                if verbose:
                    print(f"  source failed: '{current_source}' ({err})")
                continue

            if verbose:
                print(f"  selected source: '{current_source}'")
            # Destination errors are local failures and must not trigger fallback.
            shutil.copytree(folder_path, path, ignore_dangling_symlinks=True)
            return

    attempts = "; ".join(
        f"{failed_source}: {type(err).__name__}: {err}"
        for failed_source, err in failures
    )
    msg = (
        f"mhm-download: failed to download test domain '{domain}' at ref "
        f"'{branch}'. Attempts: {attempts}"
    )
    raise RuntimeError(msg) from failures[-1][1]


def cli(argv=None):  # pragma: no cover
    """Command line interface to download test domains for mHM."""
    from .. import __version__

    parser = argparse.ArgumentParser(
        description="Download tool to retrieve the test domains for mHM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=__version__,
        help="display version information",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="be verbose")
    parser.add_argument(
        "-b",
        "--branch",
        help=(
            "branch, tag, or commit of the mHM repository to take the "
            "test domain from, by default tag determined from the mHM version"
        ),
    )
    parser.add_argument(
        "-d",
        "--domain",
        type=int,
        default=1,
        choices=VALID_DOMAINS,
        help="test domain '1' or '2'",
    )
    parser.add_argument(
        "-p",
        "--path",
        help=(
            "destination path for the downloaded folder, "
            "by default the original folder name in the current directory"
        ),
    )
    parser.add_argument(
        "-s",
        "--source",
        default="auto",
        choices=VALID_SOURCES,
        help="download source; 'auto' tries git.ufz first and then github",
    )
    # parse arguments
    args = parser.parse_args(argv)
    download_test(
        branch=args.branch,
        domain=args.domain,
        path=args.path,
        verbose=args.verbose,
        source=args.source,
    )
