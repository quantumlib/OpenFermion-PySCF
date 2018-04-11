#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""Tests many modules to call pyscf functions."""
from __future__ import absolute_import

from openfermion.hamiltonians import MolecularData
from openfermionpyscf import run_pyscf
from openfermionpyscf import PyscfMolecularData


geometry = [('H', (0., 0., 0.)), ('H', (0., 0., 0.7414))]
basis = '6-31g'
multiplicity = 1
charge = 0
molecule = MolecularData(geometry,
                         basis,
                         multiplicity,
                         charge)


def test_run_pyscf():
    new_mole = run_pyscf(molecule,
                         run_scf=True,
                         run_mp2=True,
                         run_cisd=True,
                         run_ccsd=True,
                         run_fci=True,
                         verbose=1)
    assert isinstance(new_mole, PyscfMolecularData)
