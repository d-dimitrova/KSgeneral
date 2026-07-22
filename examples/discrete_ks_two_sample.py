"""Run a two-sample discrete KS test through the KSgeneral ctypes binding.

Build the native library first, for example:

    cmake -S . -B build
    cmake --build build
    python examples/discrete_ks_two_sample.py build/libksgeneral.so

On Windows, pass the path to ``ksgeneral.dll`` instead.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

dll_directory_handles = None
if sys.platform == "win32" and sys.version_info >= (3, 8):
    dll_directory_handles = []
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        directory = directory.strip().strip('"')

        if directory and os.path.isdir(directory):
            try:
                handle = os.add_dll_directory(directory)
                dll_directory_handles.append(handle)
            except OSError:
                pass

from python.ksgeneral_ctypes import load


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "library",
        type=Path,
        help="Path to libksgeneral.so, libksgeneral.dylib, or ksgeneral.dll",
    )
    args = parser.parse_args()

    # Integer-valued observations intentionally include ties, so this exercises
    # the discrete/tied-data preprocessing before calling the native library.
    sample_x = [0, 0, 0, 1]
    sample_y = [1, 1, 2, 2]

    ksg = load(args.library)
    result = ksg.two_sample_ks_test(sample_x, sample_y, alternative="two-sided")

    print(f"KSgeneral ABI: {ksg.abi_version} ({ksg.version})")
    print(f"D statistic: {result.statistic:.12g}")
    print(f"p-value: {result.p_value:.12g}")
    print(f"tie multiplicities: {result.multiplicities}")


if __name__ == "__main__":
    main()
