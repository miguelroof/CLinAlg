#!/usr/bin/python3

import os
import fnmatch
import sysconfig
from setuptools import setup, find_packages, Extension
from setuptools.command.build_py import build_py as _build_py

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


class build_py(_build_py):

    def find_package_modules(self, package, package_dir):
        ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')
        modules = super().find_package_modules(package, package_dir)
        filtered_modules = []
        for (pkg, mod, filepath) in modules:
            if os.path.exists(filepath.replace('.py', ext_suffix)):
                continue
            filtered_modules.append((pkg, mod, filepath,))
        return filtered_modules


# extensions = [
#     Extension("cypack.utils", ["src/cypack/utils.pyx"]),
#     Extension("cypack.answer", ["src/cypack/answer.pyx"]),
#     Extension("cypack.fibonacci", ["src/cypack/fibonacci.pyx"]),
#     Extension(
#         "cypack.sub.wrong",
#         ["src/cypack/sub/wrong.pyx", "src/cypack/sub/helper.c"]
#     ),
# ]

extensions = [
    Extension("CLinAlg.AlgTool", ["src/CLinAlg/AlgTool.pyx"]),
    Extension("CLinAlg.AlgVector", ["src/CLinAlg/AlgVector.pyx"]),
    Extension("CLinAlg.AlgMatrix", ["src/CLinAlg/AlgMatrix.pyx"]),
    Extension("CLinAlg.AlgQuaternion", ["src/CLinAlg/AlgQuaternion.pyx"]),
    Extension("CLinAlg.AlgTransformation", ["src/CLinAlg/AlgTransformation.pyx"]),
Extension("CLinAlg.AlgLine", ["src/CLinAlg/AlgLine.pyx"]),
]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0))) and cythonize is not None
CYTHONIZE = True

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)

with open("requirements.txt") as fp:
    install_requires = fp.read().strip().split("\n")

# with open("requirements-dev.txt") as fp:
#     dev_requires = fp.read().strip().split("\n")

setup(
    ext_modules=extensions,
    install_requires=install_requires,
    cmdclass={'build_py': build_py}
    # extras_require={
    #     "dev": dev_requires,
    #     "docs": ["sphinx", "sphinx-rtd-theme"]
    # },
)
