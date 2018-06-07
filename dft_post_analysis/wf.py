#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Calculate and return work function for periodic slab.

Author(s): I got the script from Colin but I think Karen wrote these methods
"""

#| - Import Modules
# from ase.units import Bohr
# from ase.io import read
import numpy as np
import os
#__|

def find_max_empty_space(atoms, edir=3):
        """Return scaled midpoint coordinate of largest empty space.

        Assuming periodic boundary conditions, finds the largest
        continuous segment of free, unoccupied space and returns its midpoint
        in scaled coordinates (0 to 1) in the edir direction (default z).
        """
        #| - find_max_empty_space
        # 0-indexed direction
        position_array = atoms.get_scaled_positions()[..., edir - 1]
        position_array.sort()
        differences = np.diff(position_array)
        # through the PBC
        differences = np.append(
            differences,
            position_array[0] + 1 - position_array[-1],
            )
        max_diff_index = np.argmax(differences)
        if max_diff_index == len(position_array) - 1:
            out = (position_array[0] + 1 + position_array[-1]) / 2 % 1
            return(out)
        else:
            out = (position_array[max_diff_index] +
                position_array[max_diff_index + 1]) / 2.
            return(out)
        #__|

def calc_wf(atoms, outdir):
    """Calculate work function of slab.

    Args:
        atoms:
        outdir:
    """
    #| - calc_wf
    hartree = 27.21138505
    rydberg = 0.5 * hartree
    bohr = 0.52917721092

    # rydberg_over_bohr = rydberg / bohr

    # Pick a good place to sample vacuum level
    edir = 3
    cell_length = atoms.cell[edir - 1][edir - 1] / bohr

    vacuum_pos_raw = find_max_empty_space(atoms, edir)

    vacuum_pos = vacuum_pos_raw * cell_length
    avg_out = open("%s/avg.out" % outdir, "r")
    record = False
    average_data = []
    lines = list(avg_out)
    for line in lines:
            # just start reading it in at 0
            if len(line.split()) == 3 and line.split()[0] == "0.000000000":
                    record = True
            elif len(line.split()) == 0:
                    record = False
            if record is True:
                    average_data.append([float(i) for i in line.split()])

    #vacuum_energy = average_data[np.abs(np.array(average_data)[..., 0] -
    # vacuum_pos).argmin()][2]

    # Get the latest Fermi energy
    fermi_data = os.popen('grep -n "Fermi" %s/log | tail -1' % outdir, "r")
    fermi_energy = float(fermi_data.readline().split()[-2])
    fermi_data.close()
    eopreg = 0.025

    # we use cell_length*eopreg*3 here since the work functions seem to
    # converge at that distance rather than *1 or *2
    vac1guess = vacuum_pos - cell_length * eopreg * 2.5
    vac2guess = vacuum_pos + cell_length * eopreg * 2.5
    if vac1guess < 0:
            #look for where we're closest to zero
            vac_pos1 = np.abs(cell_length + vacuum_pos - cell_length *
                eopreg * 2.5 - np.array(average_data)[..., 0]).argmin()
            #print 'shifted vac1 for %s'%outdir
    else:
            vac_pos1 = np.abs(np.array(average_data)[..., 0] -
                vacuum_pos + cell_length * eopreg * 2.5).argmin()
            #print 'normal vac1'
    if vac2guess > cell_length:
            vac_pos2 = np.abs(vacuum_pos + cell_length * eopreg * 2.5 -
                cell_length - np.array(average_data)[..., 0]).argmin()
            #print 'shifted vac2 for %s'%outdir
    else:
            vac_pos2 = np.abs(np.array(average_data)[..., 0] - vacuum_pos -
                cell_length * eopreg * 2.5).argmin()
            #print 'normal vac2'
    vacuum_energy1 = average_data[vac_pos1][1]
    vacuum_energy2 = average_data[vac_pos2][1]
    wf = [
        vacuum_energy1 * rydberg - fermi_energy,
        vacuum_energy2 * rydberg - fermi_energy,
        ]

    return(wf)
    #__|
