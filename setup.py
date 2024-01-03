#!/usr/bin/env python
"""
local_volume_database: database of the properties of Local Volume galaxies and star clusters.

Project website: https://github.com/apace7/local_volume_database

"""
import os

from setuptools import setup

_name = "local_volume_database"
_version = "0.0.1"

setup(
    name=_name,
    version=_version,
    author='Andrew B. Pace',
    author_email='pvpace1@gmail.com',
    packages=[],  # 'local_volume_database',
    package_data={
        'data/': ['*.csv'],
    },
    url='https://github.com/apace7/local_volume_database',
    license='LICENSE',
    description=
    'Database of the properties of Local Volume dwarf galaxies and Milky Way star cluster',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy",
        "astropy",
    ],
)
