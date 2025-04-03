import os
import numpy
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize

set_gcc_args   = (os.getenv("PYR_SET_GCC_ARGS")   == "1")
profile_cython = (os.getenv("PYR_PROFILE_CYTHON") == "1")

class custom_build_ext(build_ext):
    def build_extensions(self):
        if set_gcc_args:
            # This removes the "default" compiler flags that would
            # otherwise get passed on to to the compiler, i.e.,
            # sysconfig.get_config_var("CFLAGS").
            self.compiler.set_executable("compiler_so", "gcc")
        build_ext.build_extensions(self)

    def build_extension(self, extension):
        if set_gcc_args:
            extension.extra_compile_args = [
                  # Warnings
                  "-Wall",
                  "-Wno-maybe-uninitialized",
                  "-Wno-unused-result",
                  "-Wno-unused-function",
                  "-Wsign-compare",
                  # Optimization
                  "-O3",
                  "-fno-semantic-interposition",
                  # Debugging
                  "-g",
                  "-DNDEBUG",
                  "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION",
                  # Code Generation
                  "-fwrapv",
                  "-fPIC",
            ]
        build_ext.build_extension(self, extension)

extensions = [
    Extension("*", ["src/pyretechnics/*.py"], include_dirs=[numpy.get_include()])
]

setup(name="pyretechnics",
      cmdclass={"build_ext": custom_build_ext},
      ext_modules=cythonize(extensions,
                            exclude="src/pyretechnics/py_types.py",
                            annotate=True,
                            compiler_directives={
                                "language_level": "3",
                                "profile": profile_cython,
                                "initializedcheck": False,
                                "cdivision": True,
                                "wraparound": False,
                                "boundscheck": False,
                            }))
