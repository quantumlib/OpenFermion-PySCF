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

"""These functions compare properties of different molecules."""
import matplotlib.pyplot
import numpy
import warnings

from openfermion.chem import (make_atom, make_atomic_ring,
                              MolecularData, periodic_table)


def latex_name(molecule):
    """Write the name of the molecule in LaTeX.

    Returns:
        name: A string giving the name in LaTeX.
    """
    # Get sorted atom vector.
    atoms = [item[0] for item in molecule.geometry]
    atom_charge_info = [(atom, atoms.count(atom)) for atom in set(atoms)]
    sorted_info = sorted(atom_charge_info,
                         key=lambda atom: molecular_data.
                         _PERIODIC_HASH_TABLE[atom[0]])

    # Name molecule and return.
    name = '{}$_{}$'.format(sorted_info[0][0], sorted_info[0][1])
    for info in sorted_info[1::]:
        name += '{}$_{}$'.format(info[0], info[1])
    return name


# Run.
if __name__ == '__main__':

    # Set plot parameters.
    matplotlib.pyplot.rcParams['text.usetex'] = True
    matplotlib.pyplot.rcParams['text.latex.unicode'] = True
    matplotlib.pyplot.rc('text', usetex=True)
    matplotlib.pyplot.rc('font', family='sans=serif')
    marker_size = 6
    line_width = 2
    axis_size = 12
    font_size = 16
    x_log = 0
    y_log = 0

    # Set chemical series parameters.
    max_electrons = 10
    spacing = 0.7414
    basis = 'sto-3g'

    # Get chemical series.
    molecular_series = []
    for n_electrons in range(2, max_electrons + 1):
        molecule = make_atomic_ring(n_electrons, spacing, basis)
        molecule.load()
        molecular_series += [molecule]

    # Get plot data.
    x_values = []
    y_values = []
    for molecule in molecular_series:

        # x-axis.
        x_label = 'Number of Electrons'
        x_values += [molecule.n_electrons]

        # y-axis.
        y_label = 'MP2 Energy'
        y_values += [molecule.mp2_energy]

        # Print.
        print('\n{} for {} = {}.'.format(x_label, molecule.name, x_values[-1]))
        print('{} for {} = {}.'.format(y_label, molecule.name, y_values[-1]))

    # Plot.
    matplotlib.pyplot.figure(0)
    matplotlib.pyplot.plot(x_values, y_values, lw=0, marker='o')

    # Set log scales.
    if y_log:
        matplotlib.pyplot.yscale('log')
    if x_log:
        matplotlib.pyplot.xscale('log')

    # Finish making the plot.
    matplotlib.pyplot.xticks(size=axis_size)
    matplotlib.pyplot.yticks(size=axis_size)
    matplotlib.pyplot.xlabel(r'%s' % x_label, fontsize=font_size)
    matplotlib.pyplot.ylabel(r'%s' % y_label, fontsize=font_size)
    matplotlib.pyplot.show()
