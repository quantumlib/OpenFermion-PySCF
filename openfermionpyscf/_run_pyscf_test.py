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

import openfermion
import openfermionpyscf


def test_load_molecular_hamiltonian():
    geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., 1.4))]

    lih_hamiltonian = openfermionpyscf.generate_molecular_hamiltonian(
            geometry, 'sto-3g', 1, 0, 2, 2)
    assert openfermion.count_qubits(lih_hamiltonian) == 4

    lih_hamiltonian = openfermionpyscf.generate_molecular_hamiltonian(
            geometry, 'sto-3g', 1, 0, 2, 3)
    assert openfermion.count_qubits(lih_hamiltonian) == 6

    lih_hamiltonian = openfermionpyscf.generate_molecular_hamiltonian(
            geometry, 'sto-3g', 1, 0, None, None)
    assert openfermion.count_qubits(lih_hamiltonian) == 12
