"""Class defining methods to extract/manipulate data in vasp raman job folders."""

#| - Import Modules
# from dft_job_automat.job_setup import DFT_Jobs_Setup
# from aws.aws_class import AWS_Queues
import pandas as pd
pd.options.display.max_colwidth = 250

import os
# from ase import io
import numpy as np
# import pickle
# import boto3
# import subprocess
#__|

class Raman_Vasp():
    """Summary line.
    TEMP
    """

    def __init__(self, methods_to_run=[]):
        """TMP_docstring.
        TEMP TEMP

        Args:
        """
        #| - __init__
        self.tmp = 42
        self.methods_to_run = methods_to_run
        #__|

    def tmp_meth(self, path_i):
        """
        """
        #| - tmp_meth
        return("tmp - tmp_meth")
        #__|

    def read_modes_file(self, path_i):
        """
        """
        #| - read_modes_file
        # raman_dat_file = path_i + "/simulation/vasp_raman.dat"
        # print("3:" + path_i)

        raman_dat_file = path_i + "/vasp_raman.dat"

        # print(raman_dat_file)

        data_list = []
        if os.path.isfile(raman_dat_file):
            # print("*(^)")
            try:
                with open(raman_dat_file) as fle:
                    lines = fle.read().splitlines()

                    if not lines:
                        data_list = np.nan
                    elif len(lines) == 1:
                        data_list = np.nan

                    else:
                        # print(len(lines))
                        for line_num, line in enumerate(lines):
                            if line_num == 0:
                                continue

                            line = line.split()

                            line_data_dict = {}
                            line_data_dict["mode_number"]   = int(line[0])
                            line_data_dict["frequency"]     = float(line[1])
                            line_data_dict["alpha"]         = float(line[2])
                            line_data_dict["beta2"]         = float(line[3])
                            line_data_dict["activity"]      = float(line[4])

                            data_list.append(line_data_dict)

            except:
                # print("##@!@!@")
                data_list = np.nan
                print("Couldn't open vasp_raman.dat")

        else:
            # print("&&*&(())")
            data_list = np.nan

        # Can't handle multiple modes at the same time right now.
        if pd.isnull(data_list) == False:
            assert len(data_list) == 1

            if len(data_list) == 1:
                data_list = data_list[0]

        return(data_list)

        #__|
