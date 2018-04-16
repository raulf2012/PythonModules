"""Class defining methods to extract/manipulate data in vasp raman job folders.

Development Notes:
    TODO Figure out how to pass additional parameters to these methods
"""

#| - Import Modules
import os

from ase import io
from ase.io.trajectory import Trajectory

import pandas as pd
import numpy as np

# My Modules
from ase_modules.ase_methods import number_of_atoms
#__|

class DFT_Methods():
    """Summary line.
    TEMP
    """
    #| - DFT_Methods **********************************************************
    def __init__(self, methods_to_run=[]):
        """Initialize DFT_Methods instance with methods_to_run list.

        Args:
            methods_to_run:
        """
        #| - __init__
        self.methods_to_run = methods_to_run
        #__|

    def gibbs_energy(self, path_i):
        """Read gibbs free energy from file.

        Included raw electronic energy and whatever correctios were applied.

        Args:
            path_i
        """
        #| - gibbs_energy
        gibbs_e = None

        fle_name = "g_energy.out"
        if os.path.exists(path_i + "/" + fle_name):
            with open(path_i + "/" + fle_name, "r") as fle:
                gibbs_e = float(fle.read().strip())

        return(gibbs_e)
        #__|

    def elec_energy(self, path_i, atoms_file="out_opt.traj"):
        """Read electronic energy from ASE atoms object.

        Args:
            path_i:
            atoms_file:
        """
        #| - elec_energy

        #| - TEMP | VASP Run Where I Write Energy and Force TO File | 180120
        # # path_i = self.root_dir + "/" + path + "/simulation/elec_energies.out"
        #
        # path_i = path
        #
        # with open(path_i) as  f:
        #     content = f.readlines()
        #
        # content = [x.strip() for x in content]
        #
        # energy = float(content[1].split(" ")[0])
        # force = float(content[1].split(" ")[1])
        #
        # return(energy)
        #__|

        # atoms_file_names = ["out_opt.traj", "out.traj"]
        #
        # for file_name in atoms_file_names:
        #     try:
        #         atoms = io.read(path + "/" + file_name)
        #         break
        #     except:
        #         pass

        atoms = self.atoms_object(path_i)[-1]
        energy = atoms.get_potential_energy()

        return(energy)
        #__|

    def atoms_object(self, path_i):
        """

        Args:
            path_i:
        """
        #| - atoms_object
        atoms_file_names = ["out_opt.traj", "out.traj"]
        for file_name in atoms_file_names:
            try:
                # traj = io.read(path_i + "/" + file_name)
                traj = Trajectory(path_i + "/" + file_name)
                # Trajectory object can't be pickled, this is a work around for now
                traj_list = []
                for image in traj:
                    traj_list.append(image)
                traj = traj_list
                break

            except:
                traj = None
                # pass

        return(traj)
        #__|

    def init_atoms(self, path_i):
        """
        Args:
            path_i:
        """
        #| - init_atoms
        atoms_file_names = ["init.traj", "init.POSCAR"]
        for file_name in atoms_file_names:
            try:
                traj = io.read(path_i + "/" + file_name)
                break

            except:
                pass

        return(traj)
        #__|

    def parse_error_file(self, path_i):
        """Parse QE ase-espresso error file for keywords indicating failure.

        TODO Move this from job_analysis to here


        Args:
            path_i:
        """
        #| - parse_error_file
        tmp = 42


        #__|

    def atom_type_num_dict(self, path_i):
        """

        Args:
            path_i:
        """
        #| - number_of_atoms
        atoms = None
        atoms_file_names = ["out_opt.traj", "out.traj"]
        for file_name in atoms_file_names:
            try:
                # print(path_i + "/" + file_name)
                atoms = io.read(path_i + "/" + file_name)
                break

            except:
                pass

        # print(atoms)
        # atoms = io.read(path_i + "/out.traj")

        out_dict = number_of_atoms(atoms)

        return([out_dict])
        #__|

    #__| **********************************************************************
