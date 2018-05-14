#!/usr/bin/env python

"""Plot projected density of states data.

Author: Raul A. Flores

Development Notes:
    TODO Include the atoms object so that I can reference atom types
"""

#| - Import Modules
# import sys
# import pickle as pickle
# import numpy as np
import pandas as pd
# from ase import io
# import plotly
import plotly.graph_objs as go
#__|

#| - Script Inputs
# filter_dict = {
#     "type": ["sum-up", "sum-down"],
#     "atom_ind": [6, 7, 8, 9, 10, 11, 12, 13, 14],
#     "band": ["d", "p"],
#     # "element": ["N"],
#     }
#
# # pdos_file = "dir_pdos/dos.pickle"
# pdos_file = "./__test__/dos.pickle"
# atoms_file = "./out_opt.traj"
#__|

#| - Methods
def plot_dos_series(
    x_data,
    y_data,
    name,
    group=None,
    ):
    """Plot DOS data series.

    Args:

    """
    #| - plot_dos_series
    trace = go.Scatter(
        x=x_data,
        y=y_data,
        name=name,
        text=name,
        fill="tozerox",
        hoverinfo="y+text",
        legendgroup=group,
        )

    return(trace)
    #__|

#__|

def plot_pdos_dos(
    pdos_data,
    filter_dict,
    atoms,
    ):
    """

    Args:
        pdos_data:
        filter_dict:
        atoms:
    """
    #| - plot_pdos_dos

    #| - Read PDOS and Atoms From File
    # try:
    #     with open(pdos_file, "r") as fle:
    #         energies, dos, pdos = pickle.load(fle)
    # except:
    #     print("No Density of States DATA Found.")
    #     sys.exit(1)

    energies, dos, pdos = pdos_data
    # atoms = io.read(atoms_file)
    #__|

    #| - Determing Whether Calclation Is Spin Polarized
    if len(dos) != 2:
        spinpol = False
    elif len(dos) == 2:
        spinpol = True
    #__|

    #| - Data Processing

    #| - Total Density of State
    # if spinpol:
    #     dos_tot_u = dos[0]
    #     dos_tot_d = dos[1]
    #     assert len(energies) == len(dos_tot_u)
    #     assert len(energies) == len(dos_tot_d)
    #
    #     # Plot Total DOS Spin Up
    #     trace = go.Scatter(
    #         x = energies,
    #         y = dos_tot_u,
    #         name = "DOS (spin up)",
    #         )
    #     data.append(trace)
    #
    #     # Plot Total DOS Spin Down
    #     trace = go.Scatter(
    #         x = energies,
    #         y = dos_tot_u,
    #         name = "DOS (spin down)",
    #         )
    #     data.append(trace)
    #
    # else:
    #     dos_tot = dos
    #     assert len(energies) == len(dos_tot)
    #
    #     # TODO Create plot trace here
    #__|

    #| - Atomic Projected Density of State
    pdos_master_data = []

    pds = plot_dos_series  # Shorten function name

    # print(spinpol)
    if spinpol:
        tmp = 42
        #| - Spinpol: True
        # for pdos_i, atom_i in zip(pdos, atoms):
        #     elem_i = atom_i.symbol
        #     ind_i = atom_i.index
        #
        #     #| - Data Format Type Dict
        #     type_dict = {}
        #
        #     type_dict["p"] = [
        #         "sum-up",
        #         "sum-down",
        #
        #         "pz-up",
        #         "pz-down",
        #
        #         "px-up",
        #         "px-down",
        #
        #         "py-up",
        #         "py-down",
        #         ]
        #
        #     type_dict["s"] = [
        #         "sum-up",
        #         "sum-down",
        #
        #         "s-up",
        #         "s-down",
        #         ]
        #
        #     type_dict["d"] = [
        #         "sum-up",
        #         "sum-down",
        #
        #         "dz2-up",
        #         "dz2-down",
        #
        #         "dzx-up",
        #         "dzx-down",
        #
        #         "dzy-up",
        #         "dzy-down",
        #
        #         "dx2-y2 up",
        #         "dx2-y2 down",
        #
        #         "dxy up",
        #         "dxy down",
        #         ]
        #
        #     #__|
        #
        #     for band_j, dos_j in pdos_i.iteritems():
        #
        #         for k_ind, type_k in enumerate(type_dict[band_j]):
        #
        #             row_i = {
        #                 "band": band_j,
        #                 "element": elem_i,
        #                 "atom_ind": ind_i,
        #                 "type": type_k,
        #                 "pdos": dos_j[k_ind],
        #                 }
        #
        #             pdos_master_data.append(row_i)
        #__|

    else:
        #| - Spinpol: False
        # TODO Plot non-spin polarized calculation
        # pass
        #| - Data Format Type Dict
        type_dict = {}
        type_dict["p"] = [
            "sum",
            "pz",
            "px",
            "py",
            ]

        type_dict["s"] = [
            "sum",
            "s",
            ]

        type_dict["d"] = [
            "sum",
            "dz2",
            "dzx",
            "dzy",
            "dx2-y2",
            "dxy",
            ]
        #__|

        pdos_master_data = []
        for pdos_i, atom_i in zip(pdos, atoms):

            elem_i = atom_i.symbol
            ind_i = atom_i.index

            for band_j, dos_j in pdos_i.items():
                tmp = 42
                for k_ind, type_k in enumerate(type_dict[band_j]):
                    row_i = {
                        "band": band_j,
                        "element": elem_i,
                        "atom_ind": ind_i,
                        "type": type_k,
                        "pdos": dos_j[k_ind],
                        }

                    pdos_master_data.append(row_i)
        #__|

    #__|

    df = pd.DataFrame(pdos_master_data)

    #__|

    # return(df)

    #| - Data Analysis
    df["name"] = df["element"] + df["atom_ind"].astype(str) + " | " + \
        df["band"] + df["type"]

    df["max_dens"] = df["pdos"].map(lambda x: x.max())
    max_dens = df["max_dens"].max()

    # print(list(df["type"]))

    # Filter Data
    for key, value in filter_dict.items():
        df = df[df[key].isin(value)]

    #__|

    #| - Plotting

    #| - Plot Settings
    plot_title_size = 20
    tick_lab_size = 16
    axes_lab_size = 18
    legend_size = 18
    #__|

    #| - Plot Layout
    layout = {
        "title": "PDOS of Fe-supported, N-doped graphene",
        "font": {
            "family": "Courier New, monospace",
            "size": plot_title_size,
            "color": "black",
            },

        #| - Axes -----------------------------------------------------------------
        "yaxis": {
            "title": "E - E<sub>fermi</sub> [eV]",
            "zeroline": True,
            "titlefont": dict(size=axes_lab_size),
            "showgrid": False,
            "range": [-6, 3],

            "tickfont": dict(
                size=tick_lab_size,
                ),
            # gridwidth = 2,
            },

        "xaxis": {
            "title": "Density of States",
            "zeroline": False,
            "titlefont": dict(size=axes_lab_size),

            "tickfont": dict(
                size=tick_lab_size,
                ),
            "showticklabels": False,
            "range": [0, max_dens],
            },
        #__| ----------------------------------------------------------------------

        #| - Legend ---------------------------------------------------------------
        "legend": {
            "traceorder": "normal",
            "font": dict(size=legend_size)
            },

        #__| ----------------------------------------------------------------------

        #| - Plot Size
        # "width": 200 * 4.,
        # "height": 200 * 3.,
        #__|

        }
    #__|

    #| - Plotting Data Series
    data = []
    for index, row in df.iterrows():
        data_i = plot_dos_series(
            row["pdos"],
            energies,
            row["name"],
            )

        data.append(data_i)
    #__|

    #__|

    return(data, layout)

    # plotly.offline.plot(
    #     {
    #         "data": data,
    #         "layout": layout,
    #         },
    #     filename="pl_pdos.html"
    #     )

    #__|


#| - __old__

# def plot_pdos_dos(
#     pdos_data,
#     filter_dict,
#     atoms,
#     ):
#     """
#
#     Args:
#         pdos_data:
#         filter_dict:
#         atoms:
#     """
#     #| - plot_pdos_dos
#     tmp = 42
#     # #| - Import Modules
#     # import sys
#     # import pickle as pickle
#     # import numpy as np
#     # import pandas as pd
#     # from ase import io
#     #
#     # import plotly
#     # import plotly.graph_objs as go
#     # #__|
#     #
#     # #| - Script Inputs
#     # # filter_dict = {
#     # #     "type": ["sum-up", "sum-down"],
#     # #     "atom_ind": [6, 7, 8, 9, 10, 11, 12, 13, 14],
#     # #     "band": ["d", "p"],
#     # #     # "element": ["N"],
#     # #     }
#     #
#     # # pdos_file = "dir_pdos/dos.pickle"
#     # pdos_file = "./__test__/dos.pickle"
#     # atoms_file = "./out_opt.traj"
#     # #__|
#     #
#     # #| - Methods
#     # def plot_dos_series(
#     #     x_data,
#     #     y_data,
#     #     name,
#     #     group=None,
#     #     ):
#     #     """Plot DOS data series.
#     #
#     #     Args:
#     #
#     #     """
#     #     #| - plot_dos_series
#     #
#     #     # trace = go.Scatter(
#     #     trace = go.Scattergl(
#     #         x=x_data,
#     #         y=y_data,
#     #         name=name,
#     #         text=name,
#     #         fill="tozerox",
#     #         hoverinfo="y+text",
#     #         legendgroup=group,
#     #         )
#     #
#     #     return(trace)
#     #     #__|
#     #
#     # #__|
#     #__|

#__|
