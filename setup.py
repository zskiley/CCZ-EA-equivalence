from __future__ import annotations

import os
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class BuildCCZBindings(build_ext):
    def run(self) -> None:
        repo_root = Path(__file__).resolve().parent
        build_script = repo_root / "python" / "build.py"
        if not build_script.is_file():
            raise RuntimeError(f"Missing build script: {build_script}")

        subprocess.run([sys.executable, str(build_script)], cwd=repo_root, check=True)

        ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
        if not ext_suffix:
            ext_suffix = ".pyd" if os.name == "nt" else ".so"

        py_dir = repo_root / "python"
        built = py_dir / f"ccz_bindings{ext_suffix}"
        if not built.is_file():
            candidates = list(py_dir.glob("ccz_bindings*.pyd")) + list(
                py_dir.glob("ccz_bindings*.so")
            )
            if not candidates:
                raise RuntimeError("Build succeeded but ccz_bindings artifact was not found")
            built = candidates[0]

        target = Path(self.get_ext_fullpath("ccz_bindings"))
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built, target)


setup(
    name="ccz",
    version="0.1.0",
    description="CCZ/EA automorphism and equivalence algorithms",
    py_modules=["auto_worker", "ccz"],
    package_dir={"": "python"},
    ext_modules=[Extension("ccz_bindings", sources=[])],
    cmdclass={"build_ext": BuildCCZBindings},
    zip_safe=False,
)
