#| - IMPORT MODULES
import os
import sys

from ase import io
import shutil
#__|

class DataFrame_Methods():
    """Summary line.
    """
    #| - DataFrame_Methods ****************************************************
    def __init__(self, dataframe):
        """
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

        Args:
            outdir:
                Places atoms objects in folder named 'atoms_objects'
        """
        #| - create_atoms_objects
        df = self.df

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        for index, row in df.iterrows():
            atoms_i = row[atoms_row]
            if type(atoms_i) == list:
                atoms_i = atoms_i[image]

            if atoms_i != None:
                file_name = row["path"][5:].replace("/", "_") + ".traj"
                atoms_i.write(file_name)
                print(file_name)
                shutil.move(file_name, outdir + "/" + file_name)

            else:
                pass
        #__|


    #__| **********************************************************************
