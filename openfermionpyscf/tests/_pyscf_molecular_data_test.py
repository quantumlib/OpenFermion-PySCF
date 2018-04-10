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

"""Tests for molecular_data."""
from __future__ import absolute_import

import numpy

from openfermionpyscf import prepare_pyscf_molecule
from openfermionpyscf import PyscfMolecularData
from pyscf import scf, mp, ci, cc, fci


geometry = [('H', (0., 0., 0.)), ('H', (0., 0., 0.7414))]
basis = '6-31g'
multiplicity = 1
charge = 0
molecule = PyscfMolecularData(geometry,
                              basis,
                              multiplicity,
                              charge)

mol = prepare_pyscf_molecule(molecule)
mol.verbose = 0
molecule._pyscf_data['mol'] = mol
molecule._pyscf_data['scf'] = mf = scf.RHF(mol).run()
molecule._pyscf_data['mp2'] = mp.MP2(mf).run()
molecule._pyscf_data['cisd'] = ci.CISD(mf).run()
molecule._pyscf_data['ccsd'] = cc.CCSD(mf).run()
molecule._pyscf_data['fci'] = fci.FCI(mf).run()


def test_accessing_rdm():
    mo = molecule.canonical_orbitals
    overlap = molecule.overlap_integrals
    h1 = molecule.one_body_integrals
    h2 = molecule.two_body_integrals
    mf = molecule._pyscf_data['scf']
    e_core = mf.energy_nuc()

    rdm1 = molecule.cisd_one_rdm
    rdm2 = molecule.cisd_two_rdm
    e_ref = molecule._pyscf_data['cisd'].e_tot
    e_tot = (numpy.einsum('pq,pq', h1, rdm1) +
             numpy.einsum('pqrs,pqrs', h2, rdm2) * .5 + e_core)
    numpy.testing.assert_almost_equal(e_tot, e_ref, 9)

    rdm1 = molecule.ccsd_one_rdm
    rdm2 = molecule.ccsd_two_rdm
    e_ref = molecule._pyscf_data['ccsd'].e_tot
    e_tot = (numpy.einsum('pq,pq', h1, rdm1) +
             numpy.einsum('pqrs,pqrs', h2, rdm2) * .5 + e_core)
    numpy.testing.assert_almost_equal(e_tot, e_ref, 7)

    rdm1 = molecule.fci_one_rdm
    rdm2 = molecule.fci_two_rdm
    #e_ref = molecule._pyscf_data['fci'].e_tot
    e_tot = (numpy.einsum('pq,pq', h1, rdm1) +
             numpy.einsum('pqrs,pqrs', h2, rdm2) * .5 + e_core)
    numpy.testing.assert_almost_equal(e_tot, -1.1516827321, 9)

    ccsd_t1 = molecule.ccsd_single_amps
    ccsd_t2 = molecule.ccsd_double_amps


