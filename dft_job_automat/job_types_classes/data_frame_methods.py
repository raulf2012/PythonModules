"""Dataframe methods.

Development Notes:

"""

#| - IMPORT MODULES
import os
import shutil

# import sys
# from ase import io
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
        atoms_row="init_atoms",
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

        for index, row in df.iterrows():
            atoms_i = row[atoms_row]
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
