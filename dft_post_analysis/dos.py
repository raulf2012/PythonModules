#!/usr/bin/env python

"""Plot projected density of states data.

Author: Raul A. Flores

Development Notes:
    TODO Include the atoms object so that I can reference atom types
"""

#| - Import Modules
import copy
from itertools import compress

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from misc_modules.numpy_methods import make_filter_list
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

def filter_pdos_data(pdos_data, percent_keep=0.4):
    """Filter dos and pdos data series to lower memory cost.

    Args:
        pdos_data:
        percent_keep:
            Fraction of data to keep, the rest is discarded
    """
    #| - filter_pdos_data
    len_data = len(pdos_data[0])
    filter_list = make_filter_list(len_data, percent_keep)

    # --------------------------------------------------------------------------------------

    new_pdos_data = ()

    # **************************
    new_pdos_data += (np.array(list(compress(pdos_data[0], filter_list))),)

    # **************************
    len_data = len(pdos_data[1])
    if len_data == 2:
        tuple_2 = [
            np.array(list(compress(pdos_data[1][0], filter_list))),
            np.array(list(compress(pdos_data[1][1], filter_list)))
            ]

    else:
        tuple_2 = np.array(list(compress(pdos_data[1], filter_list)))


    new_pdos_data += (tuple_2,)

    # **************************
    new_list = []
    for i_ind, atom_i in enumerate(pdos_data[2]):
        dict_i = {}
        for key, value in atom_i.items():

            series_list_new = []
            for j_ind, data_series_i in enumerate(value):
                tmp = np.array(list(compress(data_series_i, filter_list)))
                series_list_new.append(tmp)

            dict_i[key] = series_list_new

        new_list.append(dict_i)

    new_pdos_data += (new_list,)


    return(new_pdos_data)
    #__|

def plot_pdos_dos(
    pdos_data,
    atoms,
    filter_dict=None,
    group=None,
    e_range=[-6, 3],  # COMBAK
    plot_title="Projected Density of States",
    ):
    """Create PDOS plot data.

    Args:
        pdos_data:
        filter_dict:
        atoms:
    """
    #| - plot_pdos_dos
    energies, dos, pdos = pdos_data

    #| - Determing Whether Calclation Is Spin Polarized
    if len(dos) != 2:
        spinpol = False
    elif len(dos) == 2:
        spinpol = True
    #__|

    #| - Data Processing

    #| - Total Density of State
    dos_data = []
    if spinpol:
        dos_tot_u = dos[0]
        dos_tot_d = dos[1]
        assert len(energies) == len(dos_tot_u)
        assert len(energies) == len(dos_tot_d)

        # Plot Total DOS Spin Up
        trace = go.Scatter(
            x=dos_tot_u,
            y=energies,
            name="DOS (spin up)",
            fill="tozerox",
            hoverinfo="y+text",
            )
        dos_data.append(trace)

        # Plot Total DOS Spin Down
        trace = go.Scatter(
            x=dos_tot_d,
            y=energies,
            name="DOS (spin down)",
            fill="tozerox",
            hoverinfo="y+text",
            )
        dos_data.append(trace)

    else:
        dos_tot = dos
        assert len(energies) == len(dos_tot)

        trace = go.Scatter(
            x=dos_tot,
            y=energies,

            name="DOS",
            )

        dos_data.append(trace)

    #__|

    #| - Atomic Projected Density of State
    pdos_master_data = []

    if spinpol:
        #| - Spinpol: True
        for pdos_i, atom_i in zip(pdos, atoms):
            elem_i = atom_i.symbol
            ind_i = atom_i.index

            #| - Data Format Type Dict
            type_dict = {}

            type_dict["p"] = [
                "sum-up",
                "sum-down",

                "pz-up",
                "pz-down",

                "px-up",
                "px-down",

                "py-up",
                "py-down",
                ]

            type_dict["s"] = [
                "sum-up",
                "sum-down",

                "s-up",
                "s-down",
                ]

            type_dict["d"] = [
                "sum-up",
                "sum-down",

                "dz2-up",
                "dz2-down",

                "dzx-up",
                "dzx-down",

                "dzy-up",
                "dzy-down",

                "dx2-y2 up",
                "dx2-y2 down",

                "dxy up",
                "dxy down",
                ]
            #__|

            # for band_j, dos_j in pdos_i.iteritems():
            for band_j, dos_j in pdos_i.items():

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

    else:
        #| - Spinpol: False

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

    #| - Data Analysis
    df["name"] = df["element"] + df["atom_ind"].astype(str) + " | " + \
        df["band"] + df["type"]

    df["max_dens"] = df["pdos"].map(lambda x: x.max())
    max_dens = df["max_dens"].max()

    # Filter Data
    if filter_dict is None:
        filter_dict = {}

    for key, value in filter_dict.items():
        df = df[df[key].isin(value)]

    #__|

    #| - Plotly Scatter Plot Creation
    data = []
    for index, row in df.iterrows():
        if group is not None:
            group_col = row[group]
        else:
            group_col = None

        data_i = plot_dos_series(
            row["pdos"],
            energies,
            row["name"],
            group=group_col,
            )

        data.append(data_i)

    pdos_data_out = data
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
        "title": plot_title,
        "font": {
            "family": "Courier New, monospace",
            "size": plot_title_size,
            "color": "black",
            },

        #| - Axes -------------------------------------------------------------
        "yaxis": {
            "title": "E - E<sub>fermi</sub> [eV]",
            "zeroline": True,
            "titlefont": dict(size=axes_lab_size),
            "showgrid": False,
            "range": e_range,

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
        #__| ------------------------------------------------------------------

        #| - Legend -----------------------------------------------------------
        "legend": {
            "traceorder": "normal",
            "font": dict(size=legend_size)
            },

        #__| ------------------------------------------------------------------

        #| - Plot Size
        # "width": 200 * 4.,
        # "height": 200 * 3.,
        #__|

        }
    #__|

    #__|

    return(dos_data, pdos_data_out, layout)
    #__|
