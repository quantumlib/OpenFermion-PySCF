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

from openfermion.hamiltonians import MolecularData

from openfermionpyscf import run_pyscf


if __name__ == '__main__':

    # Set chemical parameters.
    element_names = ['H', 'H']
    basis = 'sto-3g'
    charge = 0
    multiplicity = 1

    # Single point at equilibrium for testing
    spacings = [0.7414]

    # Add points for a full dissociation curve from 0.1 to 3.0 angstroms
    spacings += [0.1 * r for r in range(1, 31)]

    # Set run options
    run_scf = 1
    run_mp2 = 1
    run_cisd = 1
    run_ccsd = 1
    run_fci = 1
    verbose = 1

    # Run Diatomic Curve
    for spacing in spacings:
        description = "{}".format(spacing)
        geometry = [[element_names[0], [0, 0, 0]],
                    [element_names[1], [0, 0, spacing]]]
        molecule = MolecularData(geometry,
                                 basis,
                                 multiplicity,
                                 charge,
                                 description)

        molecule = run_pyscf(molecule,
                             run_scf=run_scf,
                             run_mp2=run_mp2,
                             run_cisd=run_cisd,
                             run_ccsd=run_ccsd,
                             run_fci=run_fci,
                             verbose=verbose)
        molecule.save()

    # Run Li H single point
    description = "1.45"
    geometry = [['Li', [0, 0, 0]],
                ['H', [0, 0, 1.45]]]
    molecule = MolecularData(geometry,
                             basis,
                             multiplicity,
                             charge,
                             description)

    molecule = run_pyscf(molecule,
                         run_scf=run_scf,
                         run_mp2=run_mp2,
                         run_cisd=run_cisd,
                         run_ccsd=run_ccsd,
                         run_fci=run_fci,
                         verbose=verbose)
    molecule.save()
