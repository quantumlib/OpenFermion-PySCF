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

"""This is a simple script for generating data."""
import os

from openfermion.hamiltonians import make_atomic_ring

from openfermionpyscf import run_pyscf


if __name__ == '__main__':

    # Set chemical parameters.
    basis = 'sto-3g'
    max_electrons = 10
    spacing = 0.7414

    # Select calculations.
    force_recompute = 1
    run_scf = 1
    run_mp2 = 1
    run_cisd = 1
    run_ccsd = 1
    run_fci = 1
    verbose = 1

    # Generate data.
    for n_electrons in range(2, max_electrons + 1):

        # Initialize.
        molecule = make_atomic_ring(n_electrons, spacing, basis)
        if os.path.exists(molecule.filename + '.hdf5'):
            molecule.load()

        # To run or not to run.
        if run_scf and not molecule.hf_energy:
            run_job = 1
        elif run_mp2 and not molecule.mp2_energy:
            run_job = 1
        elif run_cisd and not molecule.cisd_energy:
            run_job = 1
        elif run_ccsd and not molecule.ccsd_energy:
            run_job = 1
        elif run_fci and not molecule.fci_energy:
            run_job = 1
        else:
            run_job = force_recompute

        # Run.
        if run_job:
            molecule = run_pyscf(molecule,
                                 run_scf=run_scf,
                                 run_mp2=run_mp2,
                                 run_cisd=run_cisd,
                                 run_ccsd=run_ccsd,
                                 run_fci=run_fci,
                                 verbose=verbose)
            molecule.save()
