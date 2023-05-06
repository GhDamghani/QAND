#!/usr/bin/env python

from setuptools import setup, find_packages


# This setup is suitable for "python setup.py develop".

setup(name='qand',
      version='0.1',
      description='Quantitative and Algebraic Nonlinear Dynamics',
      author='Hossein Ghasem Damghani',
      license='GPL-3.0',
      zip_safe=False,
      install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "tqdm",
        "julia",
        "diffeqpy"
    ],
      packages=find_packages(),
      )
