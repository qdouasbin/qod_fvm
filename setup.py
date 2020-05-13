# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='qod_fvm',
    version='0.0.1',
    description='Quasi 1D finite volume solver for the PSAAP3 project',
    long_description=readme,
    author='Quentin Douasbin',
    author_email='quentin.douasbin@gmail.com',
    url='https://github.com/qdouasbin/1D_FV_Solver.git',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
