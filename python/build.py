#!/usr/bin/env python3
import os
import subprocess
import sys
import sysconfig
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    py_dir = Path(__file__).resolve().parent

    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    if not ext_suffix:
      ext_suffix = ".pyd" if os.name == "nt" else ".so"
    output = py_dir / f"ccz_bindings{ext_suffix}"

    include_dir = sysconfig.get_paths().get("include")
    if not include_dir:
      print("Could not determine Python include directory.", file=sys.stderr)
      return 2

    cxx = os.environ.get("CXX", "g++")

    sources = [
        "python/python.cpp",
        "algorithms.cpp",
        "ccz_auto.cpp",
        "ea_auto.cpp",
        "ccz_equivalence.cpp",
        "ea_equivalence.cpp",
        "dfs_equivalence.cpp",
        "quadratic.cpp",
        "dfs_auto.cpp",
        "partial_map.cpp",
        "partition_branch.cpp",
        "field_basics.cpp",
        "graph.cpp",
        "hyperplane.cpp",
        "ordered_partition.cpp",
        "weighted_fft.cpp",
        "refine.cpp",
        "groups/perm.cpp",
        "groups/group_ops.cpp",
        "groups/orbit_traversal.cpp",
        "groups/schreier_sims.cpp",
        "groups/graph_point_perm.cpp",
        "groups/semilinear_seed.cpp",
        "groups/orbit_candidates.cpp",
        "groups/dfs_group_helpers.cpp",
    ]
    source_paths = [str(repo_root / s) for s in sources]

    cmd = [
        cxx,
        "-std=c++17",
        "-O3",
        "-DNDEBUG",
        "-shared",
        f"-I{repo_root}",
        f"-I{include_dir}",
    ]

    if os.name != "nt":
        cmd.append("-fPIC")
    else:
        cmd.extend(["-static-libgcc", "-static-libstdc++"])

    cmd.extend(source_paths)
    if os.name == "nt":
        pyver = f"{sys.version_info.major}{sys.version_info.minor}"
        libs_dir = Path(sys.base_prefix) / "libs"
        cmd.extend([f"-L{libs_dir}", f"-lpython{pyver}"])
    cmd.extend(["-o", str(output)])

    print("Building:", output)
    print("Command:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, cwd=repo_root)
    except subprocess.CalledProcessError as exc:
        print(f"Build failed with exit code {exc.returncode}", file=sys.stderr)
        return exc.returncode

    print("Built extension:", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
