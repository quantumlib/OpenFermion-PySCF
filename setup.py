import io

from setuptools import setup, find_packages

# This reads the __version__ variable from openfermionpyscf/_version.py
exec(open('openfermionpyscf/_version.py').read())

# Readme file as long_description:
long_description = io.open('README.rst', encoding='utf-8').read()

# Read in requirements.txt
requirements = open('requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

setup(
    name='openfermionpyscf',
    version=__version__,
    author='The OpenFermion Developers',
    author_email='help@openfermion.org',
    url='http://www.openfermion.org',
    description='A plugin allowing OpenFermion to interface with PySCF.',
    long_description=long_description,
    install_requires=requirements,
    license='Apache 2',
    packages=find_packages()
)
