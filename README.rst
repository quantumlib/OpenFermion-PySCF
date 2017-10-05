OpenFermion-PySCF
=================

.. image:: https://badge.fury.io/py/openfermionpyscf.svg
    :target: https://badge.fury.io/py/openfermionpyscf

`OpenFermion <http://openfermion.org>`__ is an open source library (licensed under Apache 2) for compiling and analyzing quantum algorithms which simulate fermionic systems.
This plugin library allows the electronic structure package `PySCF <http://github.com/sunqm/pyscf>`__ (licensed under BSD-2-Clause) to interface with OpenFermion.

Installation
------------

To start using OpenFermion-Psi4, first install `OpenFermion <http://openfermion.org>`__ and
`PySCF <http://github.com/sunqm/pyscf>`__. To install the latest development version of OpenFermion-PySCF,
clone `this <http://github.com/quantumlib/OpenFermion-PySCF>`__ git repo, change directory to the top level folder and run:

.. code-block:: bash

  python -m pip install -e .

Alternatively, if using OpenFermion-PySCF as a library, one can install the last official PyPI release with:

.. code-block:: bash

  python -m pip install --pre --user openfermionpyscf

Also be sure to take a look at the ipython notebook demos in the examples folder of this repository.

Authors
-------

`Ryan Babbush <http://ryanbabbush.com>`__ (Google),
`Jarrod McClean <http://jarrodmcclean.com>`__ (Google),
`Ian Kivlichan <http://aspuru.chem.harvard.edu/ian-kivlichan/>`__ (Harvard),
Damian Steiger (ETH Zurich),
Thomas Haner (ETH Zurich),
Craig Gidney (Google) and
Dave Bacon (Google).

Disclaimer
----------
Copyright 2017 The OpenFermion Developers.
This is not an official Google product.
