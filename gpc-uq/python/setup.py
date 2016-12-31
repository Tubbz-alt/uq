""" Set up gPC iterative rotation package """
from setuptools import setup

setup(
    name='gpc',
    description="UQ library based on generalized polynomial chaos methods",
    version="0.0.1",
    packages=['gpc'],
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy'
    ],
)
