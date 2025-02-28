import os
from setuptools import setup
from Cython.Build import cythonize
import numpy

profile_cython = (os.getenv("PROFILE_CYTHON") == "1")

setup(
    name="pyretechnics",
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize(
        "src/pyretechnics/*.py", # module list
        exclude="src/pyretechnics/py_types.py",
        annotate=True,
        compiler_directives={
            "language_level": "3",
            "profile": profile_cython,
            "initializedcheck": False,
            "cdivision": True,
            "wraparound": False,
            "boundscheck": False,
        },
    ),
)
