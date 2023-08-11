#!/usr/bin/env python
"""
local_volume_database: basebase of the properties of local group galaxies and star clusters.

Project website: https://github.com/apace7/local_volume_database

"""
import os

from setuptools import setup

_name = "local_volume_database"
_version = ""

setup(
    name=_name,
    version=_version,
    author='Andrew B. Pace',
    author_email='pvpace1@gmail.com',
    packages=['local_volume_database'],
    package_data={'data/':['*.csv'],},
    url='https://github.com/apace7/local_volume_database',
    license='LICENSE',
    description='Database of the properties of local field dwarf galaxies and Milky Way Star Cluster',
    long_description=open('README.md').read(),
    install_requires=[
      "numpy",
      "astropy",
    ],
)