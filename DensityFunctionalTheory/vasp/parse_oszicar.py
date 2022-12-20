"""
"""

#| - Import  Modules
import os
import sys

import numpy as np
import pandas as pd

from pathlib import Path

# import chart_studio.plotly as py
import plotly.graph_objs as go

# ###############################################
from vasp.vasp_methods import parse_incar
#__|


def parse_oszicar(vasp_dir=".", oszicar_lines=None):
    """
    """
    #| - parse_oszicar
    root_dir = vasp_dir

    if oszicar_lines is None:
        path_i = os.path.join(root_dir, "OSZICAR")

        from pathlib import Path
        my_file = Path(path_i)
        if my_file.is_file():
            with open(path_i, "r") as f:
                oszicar_lines = f.read().splitlines()
        else:
            out_dict = dict(
                ion_step_conv_dict=None,
                N_tot=None,
                num_N_dict=None,
                )
            return(None)

    #| - Parse INCAR if available
    path_i = os.path.join(
        root_dir,
        "INCAR")
    my_file = Path(path_i)
    if my_file.is_file():
        with open(path_i, "r") as f:
            incar_lines = f.read().splitlines()

        incar_dict = parse_incar(incar_lines)

        nsw_i = incar_dict["NSW"]
        nelm_i = incar_dict["NELM"]
        incar_parsed = True
    else:
        incar_parsed = False
    #__|

    #| - First pass over OSZICAR, getting line groups
    line_beginnings = ["DAV:", "RMM:", ]

    lines_groups = []

    group_lines_i = []
    for line_i in oszicar_lines:

        if line_i[0:4] in line_beginnings:
            group_lines_i.append(line_i)

        if "F= " in line_i:
            # print("IDJIFSD")
            lines_groups.append(group_lines_i)
            group_lines_i = []

    # This should add the final group_lines in the case that it hasn't finished yet
    if len(oszicar_lines) > 0:
        if "F= " not in oszicar_lines[-1]:
            lines_groups.append(group_lines_i)
    #__|

    #| - Main loop
    N_tot = 0.

    ion_step_conv_dict = dict()
    for ion_step_i, lines_group_i in enumerate(lines_groups):

        data_dict_list = []
        for line_i in lines_group_i:
            data_dict_i = dict()

            line_list_i = [i for i in line_i.split(" ") if i != ""]

            N_i = line_list_i[1]

            if N_i == "***":
                data_dict_prev = data_dict_list[-1]
                N_prev = data_dict_prev["N"]
                N_i = N_prev + 1

            data_dict_i["N"] = int(N_i)

            E_i = line_list_i[2]
            data_dict_i["E"] = float(E_i)

            dE_i = line_list_i[3]
            data_dict_i["dE"] = float(dE_i)

            d_eps_i = line_list_i[4]
            data_dict_i["d_eps"] = float(d_eps_i)

            ncg_i = line_list_i[5]
            data_dict_i["ncg"] = int(ncg_i)

            rms_i = line_list_i[6]
            data_dict_i["rms"] = float(rms_i)

            if len(line_list_i) > 7:
                rms_c_i = line_list_i[7]
                data_dict_i["rms_c"] = float(rms_c_i)

            # #################################################
            data_dict_list.append(data_dict_i)


        if len(data_dict_list) > 0:
            df_i = pd.DataFrame(data_dict_list)
            N_tot += df_i.N.max()
            ion_step_conv_dict[ion_step_i] = df_i

    #__|


    num_N_dict = dict()
    for key, df_i in ion_step_conv_dict.items():
        num_N_i = df_i.N.max()
        num_N_dict[key] = num_N_i

    out_dict = dict(
        ion_step_conv_dict=ion_step_conv_dict,
        N_tot=N_tot,
        num_N_dict=num_N_dict,
        )
    return(out_dict)
    #__|
