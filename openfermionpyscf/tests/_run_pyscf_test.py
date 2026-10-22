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

import os
import runpy

import openfermion.chem.molecular_data
from openfermion.chem import MolecularData
import openfermion.config
from openfermionpyscf import run_pyscf
from openfermionpyscf import PyscfMolecularData


geometry = [('H', (0., 0., 0.)), ('H', (0., 0., 0.7414))]
basis = '6-31g'
multiplicity = 1
charge = 0


def test_run_pyscf(tmp_path):
    molecule = MolecularData(geometry,
                             basis,
                             multiplicity,
                             charge,
                             data_directory=str(tmp_path))
    new_mole = run_pyscf(molecule,
                         run_scf=True,
                         run_mp2=True,
                         run_cisd=True,
                         run_ccsd=True,
                         run_fci=True,
                         verbose=1)
    assert isinstance(new_mole, PyscfMolecularData)
    assert os.path.exists(f'{new_mole.filename}.hdf5')


def test_run_pyscf_without_save(tmp_path):
    test_molecule = MolecularData(geometry,
                                  basis,
                                  multiplicity,
                                  charge,
                                  data_directory=str(tmp_path))

    new_mole = run_pyscf(test_molecule,
                         run_scf=True,
                         save=False)
    assert isinstance(new_mole, PyscfMolecularData)
    assert not os.path.exists(f'{new_mole.filename}.hdf5')

    new_mole_saved = run_pyscf(test_molecule,
                               run_scf=True)
    assert isinstance(new_mole_saved, PyscfMolecularData)
    assert os.path.exists(f'{new_mole_saved.filename}.hdf5')


def test_generate_diatomic(tmp_path, monkeypatch):
    """Verify generate_diatomic.py saves files to the working directory."""
    readonly_dir = tmp_path / 'readonly_data'
    readonly_dir.mkdir()
    readonly_dir.chmod(0o555)

    work_dir = tmp_path / 'work'
    work_dir.mkdir()

    example_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', '..', 'examples',
                     'generate_diatomic.py'))

    try:
        monkeypatch.setattr(openfermion.chem.molecular_data, 'DATA_DIRECTORY',
                            str(readonly_dir))
        monkeypatch.setattr(openfermion.config, 'DATA_DIRECTORY',
                            str(readonly_dir))
        monkeypatch.chdir(work_dir)

        runpy.run_path(example_path, run_name='__main__')

        assert (work_dir / 'H2_sto-3g_singlet_0.7414.hdf5').exists()
        assert (work_dir / 'H1-Li1_sto-3g_singlet_1.45.hdf5').exists()
        assert len(list(work_dir.glob('*.hdf5'))) == 32
        assert not list(readonly_dir.iterdir())
    finally:
        readonly_dir.chmod(0o755)


def test_generate_molecular_hamiltonian_does_not_save(tmp_path, monkeypatch):
    monkeypatch.setattr(
        'openfermion.chem.molecular_data.DATA_DIRECTORY', str(tmp_path))
    geometry = [('H', (0., 0., 0.)), ('H', (0., 0., 0.7414))]
    openfermionpyscf.generate_molecular_hamiltonian(geometry, 'sto-3g', 1, 0)
    assert not any(tmp_path.iterdir())
