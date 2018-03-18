"""Class defining methods to extract/manipulate data in vasp raman job folders."""

#| - Import Modules
# from dft_job_automat.job_setup import DFT_Jobs_Setup
# from aws.aws_class import AWS_Queues

import pandas as pd

import os
from ase import io
import numpy as np

from ase.io.trajectory import Trajectory

# My Modules
from ase_modules.ase_methods import number_of_atoms

#__|

class DFT_Methods():
    """Summary line.
    TEMP
    """
    #| - DFT_Methods **********************************************************
    def __init__(self, methods_to_run=[]):
        """TMP_docstring.
        TEMP TEMP

        Args:
        """
        #| - __init__
        self.tmp = 42
        self.methods_to_run = methods_to_run
        #__|

    # def elec_energy(self, path, atoms_file="out.traj"):
    def elec_energy(self, path, atoms_file="out_opt.traj"):
        """
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

        atoms = self.atoms_object(path)[-1]
        energy = atoms.get_potential_energy()

        return(energy)
        #__|

    def atoms_object(self, path_i):
        """
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
        """
        Parse QE ase-espresso error-out file for keywords indicating failure
        """
        #| - parse_error_file
        tmp = 42



        #__|

    def atom_type_num_dict(self, path_i):
        """
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
