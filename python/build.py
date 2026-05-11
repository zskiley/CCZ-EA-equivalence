#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path


def _is_compiler_available(cxx: str) -> bool:
    return shutil.which(cxx) is not None or Path(cxx).is_file()


def _windows_compiler_candidates() -> list[str]:
    candidates = [
        "g++",
        "c++",
        "clang++",
        r"C:\msys64\ucrt64\bin\g++.exe",
        r"C:\msys64\mingw64\bin\g++.exe",
        r"C:\mingw64\bin\g++.exe",
    ]

    local_app_data = os.environ.get("LOCALAPPDATA")
    if local_app_data:
        winget_pkgs = Path(local_app_data) / "Microsoft" / "WinGet" / "Packages"
        if winget_pkgs.is_dir():
            for pkg_dir in winget_pkgs.glob("BrechtSanders.WinLibs.POSIX.UCRT*"):
                exe = pkg_dir / "mingw64" / "bin" / "g++.exe"
                if exe.is_file():
                    candidates.insert(0, str(exe))
                    break

    return candidates


def _resolve_cxx() -> str | None:
    env_cxx = os.environ.get("CXX")
    if env_cxx:
        return env_cxx if _is_compiler_available(env_cxx) else None

    candidates = (
        _windows_compiler_candidates()
        if os.name == "nt"
        else ["g++", "c++", "clang++"]
    )
    for cxx in candidates:
        if _is_compiler_available(cxx):
            return cxx
    return None


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

    requested_cxx = os.environ.get("CXX")
    cxx = _resolve_cxx()
    if cxx is None:
        if requested_cxx:
            print(
                f"Could not find C++ compiler from CXX='{requested_cxx}'",
                file=sys.stderr,
            )
        else:
            print("Could not find a usable C++ compiler.", file=sys.stderr)
        if os.name == "nt":
            print(
                "Windows fix: install MinGW g++ (for example WinLibs/MSYS2) and "
                "ensure g++.exe is on PATH, or set CXX to full g++.exe path.",
                file=sys.stderr,
            )
            print(
                "Example:\n"
                "  winget install -e --id BrechtSanders.WinLibs.POSIX.UCRT\n"
                "  set CXX=C:\\path\\to\\g++.exe",
                file=sys.stderr,
            )
        else:
            print(
                "Install a C++ compiler and ensure it is on PATH, or set CXX "
                "to its full path (for example: CXX=g++).",
                file=sys.stderr,
            )
        return 2

    sources = [
        "src/python.cpp",
        "src/algorithms.cpp",
        "src/ccz_auto.cpp",
        "src/ea_auto.cpp",
        "src/ccz_equivalence.cpp",
        "src/ea_equivalence.cpp",
        "src/dfs_equivalence.cpp",
        "src/quadratic.cpp",
        "src/dfs_auto.cpp",
        "src/ambient_affine.cpp",
        "src/partial_map.cpp",
        "src/gf2_linear.cpp",
        "src/partition_branch.cpp",
        "src/field_basics.cpp",
        "src/graph.cpp",
        "src/hyperplane.cpp",
        "src/ordered_partition.cpp",
        "src/weighted_fft.cpp",
        "src/refine.cpp",
        "src/groups/perm.cpp",
        "src/groups/group_ops.cpp",
        "src/groups/orbit_traversal.cpp",
        "src/groups/schreier_sims.cpp",
        "src/groups/graph_point_perm.cpp",
        "src/groups/orbit_candidates.cpp",
        "src/groups/dfs_group_helpers.cpp",
    ]
    source_paths = [str(repo_root / s) for s in sources]

    cmd = [
        cxx,
        "-std=c++17",
        "-O3",
        "-DNDEBUG",
        "-shared",
        f"-I{repo_root}",
        f"-I{repo_root / 'src'}",
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
    except FileNotFoundError:
        print(
            f"Failed to launch compiler '{cxx}'. Ensure it exists and is executable.",
            file=sys.stderr,
        )
        if os.name == "nt":
            print(
                "Windows fix: set CXX to full g++.exe path, then retry.",
                file=sys.stderr,
            )
        return 2
    except subprocess.CalledProcessError as exc:
        print(f"Build failed with exit code {exc.returncode}", file=sys.stderr)
        return exc.returncode

    print("Built extension:", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
