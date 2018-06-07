#!/usr/bin/env python

"""Dataframe methods.

Author: Raul A. Flores

Development Notes:
"""

#| - IMPORT MODULES
import os
import shutil
#__|

class DataFrame_Methods():
    """Summary line."""

    #| - DataFrame_Methods ****************************************************
    def __init__(self, dataframe):
        """
        Initialize with dataframe.

        Args:
            dataframe:
                Pandas dataframe object created by jobs_setup, jobs_analysis
                classes
        """
        #| - __init__
        self.df = dataframe
        #__|

    def create_atoms_objects(self,
        outdir="atoms_objects",
        # atoms_row="init_atoms",
        atoms_row="atoms_object",
        image=-1,
        ):
        """
        Create directory of atoms objects corresponding to job folders.

        Args:
            outdir:
                Places atoms objects in folder named 'atoms_objects'
            atoms_row:
            image:
        """
        #| - create_atoms_objects
        df = self.df

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        else:
            err_mess = "atoms folder already created, "
            err_mess += "delete or move folder and run command again"
            raise RuntimeError(err_mess)

        if atoms_row in list(df):
            pass
        else:
            atoms_row = "init_atoms"
            print("Couldn't find " + atoms_row + ", using 'init_atoms' instead")

        for index, row in df.iterrows():
            atoms_i = row[atoms_row]

            print(atoms_i)
            if type(atoms_i) == list:
                atoms_i = atoms_i[image]

            if atoms_i is not None:
                file_name = row["path"][5:].replace("/", "_") + ".traj"
                atoms_i.write(file_name)
                print(file_name)
                shutil.move(file_name, outdir + "/" + file_name)

            else:
                pass
        #__|

    #__| **********************************************************************
