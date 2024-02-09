"""!
Download routines to get test domains for mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
@ingroup mhm
"""

import argparse
import os
import shutil
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.error import HTTPError
from urllib.request import urlretrieve

# query = "?path={folder}"
MHM_URL = "https://git.ufz.de/mhm/mhm/-/archive/{branch}/mhm-{branch}.{format}{query}"
# when renaming develop to main, we need to be ready
BRANCH_MAP = {"main": "develop", "develop": "main"}
VALID_DOMAINS = [1, 2]


def _dl(branch, format, folder, filename):
    url = MHM_URL.format(branch=branch, format=format, query=f"?path={folder}")
    urlretrieve(url, filename)


def download_test(branch=None, domain=1, path=None, verbose=False):
    """
    Download a test domain for mHM.

    @param branch (str, optional): Branch, tag, or commit of the mHM repository to
        take the test domain from, by default tag determined from the mHM version
    @param domain (int, optional): Test domain 1 or 2, by default 1
    @param path (pathlike, optional): Destination path for the downloaded folder,
        by default original name of the test domain folder
    """
    # format fixed to tar.gz
    format = "tar.gz"
    # determine branch from mhm version
    if branch is None:
        from . import __version__

        branch = "main" if "dev" in __version__ else f"v{__version__}"
    # check test domain
    if domain not in VALID_DOMAINS:
        msg = f"mhm-download: 'domain' needs to be 1 or 2. Got: '{domain}'"
        raise ValueError(msg)
    folder = "test_domain" if domain == 1 else "test_domain_2"
    path = Path(path or folder)
    if verbose:
        print(f"Downloading mHM test domain '{domain}'")
        print(f"  branch: '{branch}'")
        print(f"  path:   '{Path(path).resolve()}'")
    # download to temporary directory
    with TemporaryDirectory() as tmp_dir:
        tar_file = Path(tmp_dir) / f"test.{format}"
        tar_dir = Path(tmp_dir) / "test"
        # redirect 'main' to 'develop' branch to be future proof
        try:
            _dl(branch, format, folder, tar_file)
        except HTTPError as err:
            if branch not in BRANCH_MAP:
                raise err
            _dl(BRANCH_MAP[branch], format, folder, tar_file)
        # extract
        with tarfile.open(tar_file, "r:gz") as tar:
            tar.extractall(path=tar_dir)
        # move sub-folder content to top level
        folder_path = tar_dir / os.listdir(tar_dir)[0] / folder
        shutil.copytree(folder_path, path, ignore_dangling_symlinks=True)


def cli(argv=None):  # pragma: no cover
    """Command line interface to download test domains for mHM."""
    from . import __version__

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
    # parse arguments
    args = parser.parse_args(argv)
    download_test(
        branch=args.branch, domain=args.domain, path=args.path, verbose=args.verbose
    )
