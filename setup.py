#!/usr/bin/python3

import os
from setuptools import setup, Extension

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


# for debuggin: Extension("CLinAlg.AlgTool", ["src/CLinAlg/AlgTool.pyx"], extra_compile_args=["-Ox", "-Zi"], extra_link_args=["-debug:full"])
extensions = [
    Extension("CLinAlg.AlgTool", ["src/CLinAlg/AlgTool.pyx"]),
    Extension("CLinAlg.AlgMatrix", ["src/CLinAlg/AlgMatrix.pyx"]),
    Extension("CLinAlg.AlgVector", ["src/CLinAlg/AlgVector.pyx"]),
    Extension("CLinAlg.AlgQuaternion", ["src/CLinAlg/AlgQuaternion.pyx"]),
    Extension("CLinAlg.AlgLine", ["src/CLinAlg/AlgLine.pyx"]),
    Extension("CLinAlg.AlgTransformation", ["src/CLinAlg/AlgTransformation.pyx"]),
    Extension("CLinAlg.AlgSegment", ["src/CLinAlg/AlgSegment.pyx"]),
    Extension("CLinAlg.AlgPlane", ["src/CLinAlg/AlgPlane.pyx"]),
    Extension("CLinAlg.AlgWire", ["src/CLinAlg/AlgWire.pyx"]),
    Extension("CLinAlg.AlgSurface", ["src/CLinAlg/AlgSurface.pyx"]),
    Extension("CLinAlg.AlgVector2D", ["src/CLinAlg/AlgVector2D.pyx"])]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 1))) and cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    # extensions = cythonize(extensions, compiler_directives=compiler_directives, gdb_debug=True)
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
    extras_require={
        "docs": ["sphinx", "sphinx-rtd-theme"]
    },
)
