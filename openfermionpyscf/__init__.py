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

"""
OpenFermion plugin to interface with PySCF.
"""

try:
    from ._pyscf_molecular_data import PyscfMolecularData
    from ._run_pyscf import prepare_pyscf_molecule, run_pyscf
except ImportError:
    raise Exception("Please install PySCF.")

from ._version import __version__
