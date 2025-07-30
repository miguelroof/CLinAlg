#!/usr/bin/python3
from setuptools import setup, find_packages
from Cython.Build import cythonize
from setuptools.extension import Extension
import sys, os

os.environ['CC'] = 'gcc'  # Asegúrate que gcc de MinGW32 está en el PATH
os.environ['LDSHARED'] = 'gcc -shared'

# Banderas de depuración
# Para GCC y Clang (Linux/macOS)
# if sys.platform != 'win32':
#     compile_args = ['-g', '-O0']  # -g para símbolos de depuración, -O0 para deshabilitar optimizaciones
#     link_args = ['-g']
# else:
#     # Para MSVC (Windows)
#     compile_args = ['/Zi', '/Od']  # /Zi para información de depuración, /Od para deshabilitar optimizaciones
#     link_args = ['/DEBUG']        # /DEBUG para el enlazador

# try:
#     from Cython.Build import cythonize
# except ImportError:
#     cythonize = None


# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
# def no_cythonize(extensions, **_ignore):
#     for extension in extensions:
#         sources = []
#         for sfile in extension.sources:
#             path, ext = os.path.splitext(sfile)
#             if ext in (".pyx", ".py"):
#                 if extension.language == "c++":
#                     ext = ".cpp"
#                 else:
#                     ext = ".c"
#                 sfile = path + ext
#             sources.append(sfile)
#         extension.sources[:] = sources
#     return extensions
# USO EL COMPILADOR DE MINGW
compile_args = ['-g', '-O0', '-m64', "-std=c11"]  # -g para símbolos de depuración, -O0 para deshabilitar optimizaciones
link_args = ['-g', '-m64']
sizeof_void_p = 8  # 64-bit
# Extension("CLinAlg", sources=["src/CLinAlg/AlgTool.pyx", "src/CLinAlg/AlgMatrix.pyx"], extra_compile_args=compile_args,
#           extra_link_args=link_args, define_macros=[('CYTHON_TRACE', '1'), ('CYTHON_TRACE_NOGIL', '1')])

# for debuggin: Extension("CLinAlg.AlgTool", ["src/CLinAlg/AlgTool.pyx"], extra_compile_args=["-Ox", "-Zi"], extra_link_args=["-debug:full"])
extensions = [

    Extension("CLinAlg.AlgTool", ["src/CLinAlg/AlgTool.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgMatrix", ["src/CLinAlg/AlgMatrix.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgVector", ["src/CLinAlg/AlgVector.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgQuaternion", ["src/CLinAlg/AlgQuaternion.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgLine", ["src/CLinAlg/AlgLine.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgTransformation", ["src/CLinAlg/AlgTransformation.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgSegment", ["src/CLinAlg/AlgSegment.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgPlane", ["src/CLinAlg/AlgPlane.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgVector2D", ["src/CLinAlg/AlgVector2D.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgWire", ["src/CLinAlg/AlgWire.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.AlgSurface", ["src/CLinAlg/AlgSurface.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.ToolFiber", ["src/CLinAlg/ToolFiber.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
    Extension("CLinAlg.ToolMath", ["src/CLinAlg/ToolMath.pyx"], extra_compile_args=compile_args,
              extra_link_args=link_args),
]

# os.putenv("Py_DEBUG", "1")
# CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 1))) and cythonize is not None
#
# if CYTHONIZE:
compiler_directives = {"language_level": 3, "embedsignature": True, 'linetrace': True, "binding": True, "profile": True}
# extensions = cythonize(extensions, compiler_directives=compiler_directives, gdb_debug=True)
extensions = cythonize(extensions, compiler_directives=compiler_directives, gdb_debug=True)
# else:
#     extensions = no_cythonize(extensions)


with open("requirements.txt") as fp:
    install_requires = fp.read().strip().split("\n")

# with open("requirements-dev.txt") as fp:
#     dev_requires = fp.read().strip().split("\n")
packages = find_packages()
setup(
    name="CLinAlg",
    version="0.5.2",
    packages=packages,
    ext_modules=extensions,
    options={'build_ext': {'compiler': 'mingw32'}},
    install_requires=install_requires,
    extras_require={
        "docs": ["sphinx", "sphinx-rtd-theme"]
    },
)
