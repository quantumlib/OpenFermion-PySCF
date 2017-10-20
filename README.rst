OpenFermion-PySCF
=================

.. image:: https://badge.fury.io/py/openfermionpyscf.svg
    :target: https://badge.fury.io/py/openfermionpyscf

.. image:: https://img.shields.io/badge/python-2.7%2C%203.4%2C%203.5%2C%203.6-brightgreen.svg

`OpenFermion <http://openfermion.org>`__ is an open source library (licensed under Apache 2) for compiling and analyzing quantum algorithms which simulate fermionic systems.
This plugin library allows the electronic structure package `PySCF <http://github.com/sunqm/pyscf>`__ (licensed under BSD-2-Clause) to interface with OpenFermion.

Installation
------------

To start using OpenFermion-PySCF, first install `PySCF
<http://github.com/sunqm/pyscf>`__.
Then, to install the latest versions of OpenFermion and OpenFermion-PySCF (in development mode):

.. code-block:: bash

  git clone https://github.com/quantumlib/OpenFermion-PySCF
  cd OpenFermion-PySCF
  python -m pip install -e .

Alternatively, to install the latest PyPI releases as libraries (in user mode):

.. code-block:: bash

  python -m pip install --user openfermionpyscf

Also be sure to take a look at the `ipython notebook demo <https://github.com/quantumlib/OpenFermion-PySCF/blob/master/examples/openfermionpyscf_demo.ipynb>`__.

How to contribute
-----------------

We'd love to accept your contributions and patches to OpenFermion-PySCF.
There are a few guidelines you need to follow.
Contributions to OpenFermion-PySCF must be accompanied by a Contributor License Agreement.
You (or your employer) retain the copyright to your contribution,
this simply gives us permission to use and redistribute your contributions as part of the project.
Head over to https://cla.developers.google.com/
to see your current agreements on file or to sign a new one.

All submissions, including submissions by project members, require review.
We use GitHub pull requests for this purpose. Consult
`GitHub Help <https://help.github.com/articles/about-pull-requests/>`__ for
more information on using pull requests.
Furthermore, please make sure your new code comes with extensive tests!
We use automatic testing to make sure all pull requests pass tests and do not
decrease overall test coverage by too much. Make sure you adhere to our style
guide. Just have a look at our code for clues. We mostly follow
`PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ and use
the corresponding `linter <https://pypi.python.org/pypi/pep8>`_ to check for it.
Code should always come with documentation.

Authors
-------

`Ryan Babbush <http://ryanbabbush.com>`__ (Google),
`Jarrod McClean <http://jarrodmcclean.com>`__ (Google),
`Ian Kivlichan <http://aspuru.chem.harvard.edu/ian-kivlichan/>`__ (Harvard),
`Damian Steiger <https://github.com/damiansteiger>`__ (ETH Zurich),
`Thomas Haener <https://github.com/thomashaener>`__ (ETH Zurich) and
`Dave Bacon <https://github.com/dabacon>`__ (Google).

How to cite
-----------
When using OpenFermion-PySCF for research projects, please cite:

* Jarrod R. McClean, Ian D. Kivlichan, Damian S. Steiger, Yudong Cao, E.
  Schuyler Fried, Craig Gidney, Thomas Häner, Vojtĕch Havlíček,
  Zhang Jiang, Matthew Neeley, Jhonathan Romero, Nicholas Rubin, Nicolas P. D.
  Sawaya, Kanav Setia, Sukin Sim, Wei Sun, Kevin Sung and Ryan Babbush.
  *OpenFermion: The Electronic Structure Package for Quantum Computers*.
  arXiv preprint. 2017.
|
* Qiming Sun, Timothy C. Berkelbach, Nick S. Blunt, George H. Booth, Sheng Guo,
  Zhendong Li, Junzi Liu, James McClain, Elvira. R. Sayfutyarova, Sandeep Sharma,
  Sebastian Wouters and Garnet Kin-Lic Chan.
  *The Python-based Simulations of Chemistry Framework (PySCF)*.
  `WIREs Compututational Molecular Science <http://onlinelibrary.wiley.com/doi/10.1002/wcms.1340/full>`__
  2017.

We are happy to include future contributors as authors on later OpenFermion releases.

Disclaimer
----------
Copyright 2017 The OpenFermion Developers.
This is not an official Google product.
