#!/usr/bin/env python

"""Class defining methods to extract/manipulate data in vasp raman job folders.

Development Notes:
    TODO Figure out how to pass additional parameters to these methods
    Make methods work for VASP and QE by using teh DFT_code attribute
"""

#| - Import Modules
import os
import pickle as pickle
from ase import io
# from ase.io.trajectory import Trajectory

# My Modules
# from ase_modules.ase_methods import number_of_atoms
from ase_modules.ase_methods import create_species_element_dict

from quantum_espresso.qe_methods import magmom_charge_data
#__|

class DFT_Methods():
    """Methods and analysis to perform within DFT jobs folders."""

    #| - DFT_Methods **********************************************************
    def __init__(self,
        methods_to_run=[],
        DFT_code="QE",  # VASP
        ):
        """Initialize DFT_Methods instance with methods_to_run list.

        Args:
            methods_to_run:
        """
        #| - __init__
        self.methods_to_run = methods_to_run
        self.DFT_code = DFT_code
        #__|

    def pdos_data(self, path_i):
        """Read pdos.pickle file and return data.

        Args:
            path_i:
        """
        #| - pdos_data
        fle_name = "dir_pdos/dos.pickle"
        if os.path.exists(path_i + "/" + fle_name):
            # with open(path_i + "/" + fle_name, "r") as fle:
            # NOTE Added "rb" & encoding="latin1" for python3 support
            with open(path_i + "/" + fle_name, "rb") as fle:
                data = pickle.load(fle, encoding="latin1")

        return(data)
        #__|

    def bands_data(self, path_i):
        """Read band_disp.pickle file and return data.

        Args:
            path_i:
        """
        #| - bands_data
        fle_name = "dir_bands/band_disp.pickle"
        # print(path_i + "/" + fle_name)
        if os.path.exists(path_i + "/" + fle_name):
            with open(path_i + "/" + fle_name, "rb") as fle:
                data = pickle.load(fle, encoding="latin1")

        return(data)
        #__|

    def magmom_charge_history(self, path_i, log="calcdir/log"):
        """Return atomic charges and magmoms thourgh SCF convergence history.

        Args:
            path_i
        """
        #| - magmom_charge_history
        df = magmom_charge_data(path_i=path_i, log=log)
        imp_col = ["atom_num", "iteration", "element"]

        magmom_history_df = df.filter(items=imp_col + ["magmom"])
        charge_history_df = df.filter(items=imp_col + ["charge"])

        out_dict = {}
        out_dict["magmom_history"] = magmom_history_df
        out_dict["charge_history"] = charge_history_df

        return(out_dict)
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

    def gibbs_correction(self, path_i):
        """Return gibbs free energy correction.

        Args:
            path_i
        """
        #| - gibbs_correction
        gibbs_corr = 0.

        fle_name = "dir_vib/gibbs_corr.out"
        if os.path.exists(path_i + "/" + fle_name):
            with open(path_i + "/" + fle_name, "r") as fle:
                gibbs_corr = float(fle.read().strip())

        return(gibbs_corr)
        #__|

    def elec_energy(self, path_i, atoms_file="out_opt.traj"):
        """Read electronic energy from ASE atoms object.

        Args:
            path_i:
            atoms_file:
        """
        #| - elec_energy
        try:
            with open(path_i + "/dir_opt/elec_e.out", "r") as fle:
                energy = float(fle.read().strip())

        except:
            pass

        # with open(path_i + "/dir_opt/elec_e.out", "r") as fle:
        #     energy = float(fle.read().strip())

        try:
            atoms = self.atoms_object(path_i)[-1]
            energy = atoms.get_potential_energy()
        except:
            pass

        try:

            atoms = io.read(path_i + "/" + "out_opt.traj")
            energy = atoms.get_potential_energy()
        except:
            pass


        return(energy)
        #__|

    def atoms_object(self, path_i):
        """Attempt to read and return atoms object.

        Args:
            path_i:
        """
        #| - atoms_object
        # atoms_file_names = ["out_opt.traj", "out.traj"]
        # 'out.traj' should be read first

        atoms_file_names = [
            "out.traj",
            "out_opt.traj",
            "OUTCAR",  # VASP
            ]

        for file_name in atoms_file_names:
            try:

                #| - try to read atoms
                if self.DFT_code == "VASP":

                    cwd = os.getcwd()
                    os.chdir(path_i)

                    traj = io.read(
                        os.path.join(path_i, file_name),
                        index=":",
                        )

                    os.chdir(cwd)

                else:

                    traj = io.read(
                        os.path.join(path_i, file_name),
                        index=":",
                        )

                break
                #__|

            except:
                traj = None

        return(traj)
        #__|


    def outcar(self, path_i):
        """Attempt to read and return outcar.

        Args:
            path_i:
        """
        #| - atoms_object
        line_list = []
        with open(os.path.join(path_i, "OUTCAR")) as fle:
            for line in fle:
                line = line.rstrip()
                line_list.append(line)

        return(line_list)

        #| - __old__
        # for file_name in atoms_file_names:
        #     try:
        #
        #         #| - try to read atoms
        #         if self.DFT_code == "VASP":
        #
        #             cwd = os.getcwd()
        #             os.chdir(path_i)
        #
        #             traj = io.read(
        #                 os.path.join(path_i, file_name),
        #                 index=":",
        #                 )
        #
        #             os.chdir(cwd)
        #
        #         else:
        #             traj = io.read(
        #                 os.path.join(path_i, file_name),
        #                 index=":",
        #                 )
        #
        #         break
        #         #__|
        #
        #     except:
        #         traj = None
        #
        # return(traj)
        #__|

        #__|


    def incar(self, path_i):
        """
        """
        #| - incar
        line_list = []
        with open(os.path.join(path_i, "INCAR")) as fle:
            for line in fle:
                line = line.rstrip()
                line_list.append(line)

        return(line_list)
        #__|

    def init_atoms(self, path_i):
        """Attempt to read and return initial atoms object.

        Args:
            path_i:
        """
        #| - init_atoms
        traj = None

        atoms_file_names = [
            "init.traj",
            "init.POSCAR",
            "out_opt.traj",
            "init.cif",
            ]

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
        print(tmp)
        #__|

    def atom_type_num_dict(self, path_i):
        """Return dictionary containing atomic count for each element.

        Args:
            path_i:
        """
        #| - atom_type_num_dict
        atoms = None
        atoms_file_names = ["out_opt.traj", "out.traj"]
        for file_name in atoms_file_names:
            try:
                print(path_i + "/" + file_name)
                atoms = io.read(path_i + "/" + file_name)
                break

            except:
                pass

        # print(atoms)
        # atoms = io.read(path_i + "/out.traj")

        # print("dft_methods")
        out_dict = create_species_element_dict(
            atoms,
            include_all_elems=False,
            elems_to_always_include=None,
            )

        print("dft_methods - atom_type_num_dict")
        print(out_dict)
        # out_dict = number_of_atoms(atoms)

        return([out_dict])
        #__|

    #__| **********************************************************************
