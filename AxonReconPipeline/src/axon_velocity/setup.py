# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

d = {}
exec(open("axon_velocity/version.py").read(), None, d)
version = d['version']
long_description = open("README.md").read()

with open("requirements.txt", mode='r') as f:
    install_requires = f.read().split('\n')
install_requires = [e for e in install_requires if len(e) > 0]

entry_points = None

setup(
    name="axon_velocity",
    version=version,
    author="Alessio Buccino, Xinyue Yuan",
    author_email="alessiop.buccino@gmail.com",
    description="Python package to reconstruct axonal branches from high-density micro-electrode array data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alejoe91/axon_velocity",
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    entry_points=entry_points,
    include_package_data=True,
)
