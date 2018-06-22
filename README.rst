OpenFermion-PySCF
=================

.. image:: https://badge.fury.io/py/openfermionpyscf.svg
    :target: https://badge.fury.io/py/openfermionpyscf

.. image:: https://travis-ci.org/quantumlib/OpenFermion-PySCF.svg?branch=master
    :target: https://travis-ci.org/quantumlib/OpenFermion-PySCF

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
`Kevin J. Sung <https://github.com/kevinsung>`__ (University of Michigan),
`Damian Steiger <https://github.com/damiansteiger>`__ (ETH Zurich),
`Dave Bacon <https://github.com/dabacon>`__ (Google),
`Yudong Cao <https://github.com/yudongcao>`__ (Harvard),
`Chengyu Dai <https://github.com/jdaaph>`__ (University of Michigan),
`E. Schuyler Fried <https://github.com/schuylerfried>`__ (Harvard),
`Craig Gidney <https://github.com/Strilanc>`__ (Google),
`Brendan Gimby <https://github.com/bgimby>`__ (University of Michigan),
`Pranav Gokhale <https://github.com/singular-value>`__ (University of Chicago),
`Thomas Häner <https://github.com/thomashaener>`__ (ETH Zurich),
`Tarini Hardikar <https://github.com/TariniHardikar>`__ (Dartmouth),
`Vojtĕch Havlíček <https://github.com/VojtaHavlicek>`__ (Oxford),
`Cupjin Huang <https://github.com/pertoX4726>`__ (University of Michigan),
`Josh Izaac <https://github.com/josh146>`__ (Xanadu),
`Zhang Jiang <https://ti.arc.nasa.gov/profile/zjiang3>`__ (NASA),
`Xinle Liu <https://github.com/sheilaliuxl>`__ (Google),
`Sam McArdle <https://github.com/sammcardle30>`__ (Oxford),
`Matthew Neeley <https://github.com/maffoo>`__ (Google),
`Thomas O'Brien <https://github.com/obriente>`__ (Leiden University),
`Isil Ozfidan <https://github.com/conta877>`__ (D-Wave Systems),
`Max Radin <https://github.com/max-radin>`__ (UC Santa Barbara),
`Jhonathan Romero <https://github.com/jromerofontalvo>`__ (Harvard),
`Nicholas Rubin <https://github.com/ncrubin>`__ (Rigetti),
`Daniel Sank <https://github.com/DanielSank>`__ (Google),
`Nicolas Sawaya <https://github.com/nicolassawaya>`__ (Harvard),
`Kanav Setia <https://github.com/kanavsetia>`__ (Dartmouth),
`Hannah Sim <https://github.com/hsim13372>`__ (Harvard),
`Mark Steudtner <https://github.com/msteudtner>`__  (Leiden University),
`Qiming Sun <https://github.com/sunqm>`__ (Caltech),
`Wei Sun <https://github.com/Spaceenter>`__ (Google),
`Daochen Wang <https://github.com/daochenw>`__ (River Lane Research),
`Chris Winkler <https://github.com/quid256>`__ (University of Chicago) and
`Fang Zhang <https://github.com/fangzh-umich>`__ (University of Michigan).

How to cite
-----------
When using OpenFermion-PySCF for research projects, please cite:

    Jarrod R. McClean, Ian D. Kivlichan, Kevin J. Sung, Damian S. Steiger,
    Yudong Cao, Chengyu Dai, E. Schuyler Fried, Craig Gidney, Brendan Gimby,
    Pranav Gokhale, Thomas Häner, Tarini Hardikar, Vojtĕch Havlíček,
    Cupjin Huang, Josh Izaac, Zhang Jiang, Xinle Liu, Matthew Neeley,
    Thomas O'Brien, Isil Ozfidan, Maxwell D. Radin, Jhonathan Romero,
    Nicholas Rubin, Nicolas P. D. Sawaya, Kanav Setia, Sukin Sim,
    Mark Steudtner, Qiming Sun, Wei Sun, Fang Zhang and Ryan Babbush.
    *OpenFermion: The Electronic Structure Package for Quantum Computers*.
    `arXiv:1710.07629 <https://arxiv.org/abs/1710.07629>`__. 2017.

as well as

    Qiming Sun, Timothy C. Berkelbach, Nick S. Blunt, George H. Booth, Sheng Guo,
    Zhendong Li, Junzi Liu, James McClain, Elvira. R. Sayfutyarova, Sandeep Sharma,
    Sebastian Wouters and Garnet Kin-Lic Chan.
    *The Python-based Simulations of Chemistry Framework (PySCF)*.
    `WIREs Compututational Molecular Science <http://onlinelibrary.wiley.com/doi/10.1002/wcms.1340/full>`__.
    2017.

We are happy to include future contributors as authors on later OpenFermion releases.

Disclaimer
----------
Copyright 2017 The OpenFermion Developers.
This is not an official Google product.
