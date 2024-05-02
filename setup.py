#!/usr/bin/env python
from setuptools import setup

setup(
    name             = "pyretechnics",
    version          = "2024.5.2",
    description      = "A Python library for simulating fire behavior in a variety of ways.",
    author           = "Gary W. Johnson",
    author_email     = "gjohnson@sig-gis.com",
    package_dir      = {"": "src"},
    install_requires = [
        "numpy",
        "rasterio",
    ],
)
