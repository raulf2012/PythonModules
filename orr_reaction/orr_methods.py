#!/usr/bin/env python

"""ORR energetics classes and methods.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import copy
import numpy as np
import pandas as pd

from plotly.graph_objs import Scatter

pd.options.mode.chained_assignment = None

from orr_reaction.orr_series import ORR_Free_E_Series
#__|

#| - __old__
def plotly_fed_layout(
    plot_title="FED",
    plot_title_size=18,
    tick_lab_size=16,
    axes_lab_size=18,
    legend_size=18,
    ):
    """
    """
    #| - plotly_fed_layout
    xax_labels = ["O2", "OOH", "O", "OH", "H2O"]
    layout = {

        "title": plot_title,

        "font": {
            "family": "Courier New, monospace",
            "size": plot_title_size,
            "color": "black",
            },

        #| - Axes --------------------------------------------------------------
        "yaxis": {
            "title": "Free Energy [eV]",
            "zeroline": True,
            "titlefont": dict(size=axes_lab_size),
            "showgrid": False,
            "tickfont": dict(
                size=tick_lab_size,
                ),
            },

        "xaxis": {
            "title": "Reaction Coordinate",
            "zeroline": True,
            "titlefont": dict(size=axes_lab_size),
            "showgrid": False,

            # "showticklabels": False,

            "ticktext": xax_labels,
            "tickvals": [1.5 * i + 0.5 for i in range(len(xax_labels))],

            "tickfont": dict(
                size=tick_lab_size,
                ),
            },
        #__| -------------------------------------------------------------------

        #| - Legend ------------------------------------------------------------
        "legend": {
            "traceorder": "normal",
            "font": dict(size=legend_size)
            },
        #__| -------------------------------------------------------------------

        #| - Plot Size
        # "width": 200 * 4.,
        # "height": 200 * 3.,
        #__|

        }

    return(layout)

    #__|

#__|


def calc_ads_e(
    df_row,
    bare_raw_e,
    correction=0.,
    oxy_ref_e=-443.70964,
    hyd_ref_e=-16.46018,
    ):
    """Calculate adsorption energies from raw DFT energetics.

    Default oxygen reference energy is based on water

    Args:
        df_row: Pandas dataframe row
        bare_raw_e: Bare slab raw DFT energy
        correction: Energy correction (ZPE, entropy, solvation, etc.)
        oxy_ref_e:
        hyd_ref_e:
    """
    #| - calc_ads_e
    row = df_row
    bare_slab = bare_raw_e
    oxy_ref = oxy_ref_e
    hyd_ref = hyd_ref_e

    #| - Oxygen & Hydrogen Atom Count

    atoms_col = "atom_type_num_dict"
    if atoms_col in list(row):
        try:
            num_O = row[atoms_col][0]["O"]
        except:
            num_O = 0

        try:
            num_H = row[atoms_col][0]["H"]
        except:
            num_H = 0

    else:

        if row["adsorbate"] == "ooh":
            num_O = 2
            num_H = 1
        elif row["adsorbate"] == "o":
            num_O = 1
            num_H = 0
        elif row["adsorbate"] == "oh":
            num_O = 1
            num_H = 1
        elif row["adsorbate"] == "bare":
            num_O = 0
            num_H = 0
    #__|

    try:
        raw_e = row["elec_energy"]
        ads_e_i = raw_e - (bare_slab + num_O * oxy_ref + num_H * hyd_ref)
        ads_e_i += correction

        # ads_e_i = raw_e - bare_slab - num_O * oxy_ref - num_H * hyd_ref
        # ads_e_i += correction
    except:
        ads_e_i = None

    return(ads_e_i)
    #__|

def df_calc_adsorption_e(
    df,

    oxy_ref,
    hyd_ref,

    bare_slab_e,
    bare_slab_var=None,

    corrections_mode="df_column",  # corr_dict
    corrections_column="gibbs_correction",

    corrections_dict=None,
    ):
    """Calculate and add adsorption energy column to data_frame.

    Args:
        df:
    """
    #| - df_calc_adsorption_e
    ads_e_list = []
    for index, row in df.iterrows():
        bare_e = bare_slab_e

        #| - Correction
        corr = 0.
        # corr = fe_corr_dict[row["adsorbate"]]

        if corrections_mode == "df_column":
            corr = row[corrections_column]

            # If "df_column" method return 0. then try to use correction_dict
            if corr == 0.:
                if corrections_dict is not None:
                    corr = corrections_dict[row["adsorbate"]]

        elif corrections_mode == "corr_dict" and corrections_dict is not None:
            corr = corrections_dict[row["adsorbate"]]
        else:
            print("No correction being applied")
            corr = 0.
        #__|

        if type(bare_slab_e) == dict:
            bare_e = bare_slab_e[row[bare_slab_var]]

        elif type(bare_slab_e) == float:
            bare_e = bare_slab_e

        ads_e_i = calc_ads_e(
            row,
            # bare_slab,
            bare_e,
            correction=corr,
            oxy_ref_e=oxy_ref,
            hyd_ref_e=hyd_ref,
            )
        ads_e_list.append(ads_e_i)

    df["ads_e"] = np.array(ads_e_list)
    #__|

def lowest_e_path(
    df,
    jobs_variables,
    color_list,
    create_ideal_series=True,
    opt_name=None,
    bias=0.,
    manual_props="*TEMP*",
    plot_title="Free Energy Diagram for the Oxygen Reduction Reaction",
    smart_format=None,
    ):
    """Find the lowest energy pathway FED.

    COMBAK

    From a set of FE pathways corresponding to different sites, the lowest
    energy states will be selected to construct a new FED.

    Args:
        df:
        jobs_variables:
            Result of Jobs.tree_level_labels
        color_list:
        bias:

    """
    #| - lowest_e_path

    #| - Grouping By Adsorbate Type
    df = copy.deepcopy(df)
    groupby = copy.deepcopy(jobs_variables)

    groupby.remove("site")
    groupby.remove("adsorbate")

    data_master = {}
    if groupby == []:
        series_list = []
        for ads_i in df.groupby("adsorbate"):

            min_e_row = ads_i[1].loc[ads_i[1]["ads_e"].idxmin()]
            series_list.append(min_e_row)

        df_i = pd.DataFrame.from_items([(s.name, s) for s in series_list]).T
        data_master[manual_props] = df_i

    else:
        for group_i in df.groupby(groupby):
            series_list = []
            for ads_i in group_i[1].groupby("adsorbate"):

                min_e_row = ads_i[1].loc[ads_i[1]["ads_e"].idxmin()]
                series_list.append(min_e_row)

            df_i = pd.DataFrame.from_items([(s.name, s) for s in series_list]).T
            data_master[group_i[0]] = df_i

    #__|

    #| - Creating Data Sets

    #| - Creating FED Datasets
    data_list = []
    # for i_cnt, (key, fe_dict) in enumerate(data_master.iteritems()):
    for i_cnt, (key, fe_dict) in enumerate(data_master.items()):

        ORR = ORR_Free_E_Plot(
            free_energy_df=fe_dict,
            )

        dat_lst = ORR.plot_fed_series(
            bias=bias,
            opt_name=opt_name,
            properties=key,
            color_list=color_list,
            i_cnt=i_cnt,
            hover_text_col="site",
            smart_format=smart_format,
            )

        data_list.extend(dat_lst)
    #__|

    #| - Creating Ideal FED Dataset
    if create_ideal_series:
        e_list_ideal = ORR.apply_bias(bias, ORR.ideal_energy)

        dat_ideal = ORR.create_plotly_series(
            e_list_ideal,
            group="Ideal",
            name="Ideal",
            color=color_list[-1],
            plot_mode="full_lines",
            )

        dat_lst = data_list + dat_ideal

    else:

        dat_lst = data_list


    #__|

    # dat_lst = data_list + dat_ideal

    #__|

    #| - Plotting

    #| - Plot Settings
    plot_title_size = 18
    tick_lab_size = 16
    axes_lab_size = 18
    legend_size = 18
    #__|

    #| - Plot Layout
    # xax_labels = ["O2", "OOH", "O", "OH", "H2O"]
    # layout = {
    #
    #     "title": plot_title,
    #
    #     "font": {
    #         "family": "Courier New, monospace",
    #         "size": plot_title_size,
    #         "color": "black",
    #         },
    #
    #     #| - Axes --------------------------------------------------------------
    #     "yaxis": {
    #         "title": "Free Energy [eV]",
    #         "zeroline": True,
    #         "titlefont": dict(size=axes_lab_size),
    #         "showgrid": False,
    #         "tickfont": dict(
    #             size=tick_lab_size,
    #             ),
    #         },
    #
    #     "xaxis": {
    #         "title": "Reaction Coordinate",
    #         "zeroline": True,
    #         "titlefont": dict(size=axes_lab_size),
    #         "showgrid": False,
    #
    #         # "showticklabels": False,
    #
    #         "ticktext": xax_labels,
    #         "tickvals": [1.5 * i + 0.5 for i in range(len(xax_labels))],
    #
    #         "tickfont": dict(
    #             size=tick_lab_size,
    #             ),
    #         },
    #     #__| -------------------------------------------------------------------
    #
    #     #| - Legend ------------------------------------------------------------
    #     "legend": {
    #         "traceorder": "normal",
    #         "font": dict(size=legend_size)
    #         },
    #     #__| -------------------------------------------------------------------
    #
    #     #| - Plot Size
    #     # "width": 200 * 4.,
    #     # "height": 200 * 3.,
    #     #__|
    #
    #     }
    #__|

    layout = plotly_fed_layout(plot_title=plot_title)

    #__|

    return(dat_lst, layout)

    #__|

def plot_all_states(
    df,
    jobs_variables,
    color_list,
    bias=0.,
    hover_text_col="site",
    create_ideal_series=True,
    plot_title="Free Energy Diagram for the Oxygen Reduction Reaction",
    smart_format=None,
    ):
    """

    Args:
        df:
        jobs_variables:
        color_list:
        bias:
        plot_title:
    """
    #| - plot_all_states

    #| - Grouping By Adsorbate Type

    groupby = copy.deepcopy(jobs_variables)
    # groupby = copy.deepcopy(Jojobs_variablesbs.tree_level_labels)
    groupby.remove("adsorbate")

    data_master = {}
    for group_i in df.groupby(groupby):

        data_master[group_i[0]] = group_i[1]
    #__|

    #| - Creating Data Sets

    #| - Creating FED Datasets
    data_list = []
    # for i_cnt, (key, fe_dict) in enumerate(data_master.iteritems()):
    for i_cnt, (key, fe_dict) in enumerate(data_master.items()):
        ORR = ORR_Free_E_Plot(
            free_energy_df=fe_dict,
            )

        dat_lst = ORR.plot_fed_series(
            bias=bias,
            properties=key,
            color_list=color_list,
            i_cnt=i_cnt,
            # hover_text_col="site"
            hover_text_col=hover_text_col,
            plot_mode="states_only",
            smart_format=smart_format,
            )

        data_list.extend(dat_lst)
    #__|

    #| - Creating Ideal FED Dataset
    if create_ideal_series:

        e_list_ideal = ORR.apply_bias(bias, ORR.ideal_energy)

        dat_ideal = ORR.create_plotly_series(
            e_list_ideal,
            group="Ideal",
            name="Ideal",
            color="red",
            plot_mode="full_lines",
            )

        dat_lst = data_list + dat_ideal
    #__|

    else:
        dat_lst = data_list

    #__|

    #| - Plotting

    #| - Plot Settings
    plot_title_size = 18
    tick_lab_size = 16
    axes_lab_size = 18
    legend_size = 12
    #__|

    #| - Plot Layout
    # xax_labels = ["O2", "OOH", "O", "OH", "H2O"]
    # layout = {
    #
    #     "title": plot_title,
    #
    #     "font": {
    #         "family": "Courier New, monospace",
    #         "size": plot_title_size,
    #         "color": "black",
    #         },
    #
    #     #| - Axes --------------------------------------------------------------
    #     "yaxis": {
    #         "title": "Free Energy [eV]",
    #         "zeroline": True,
    #         "titlefont": dict(size=axes_lab_size),
    #         "showgrid": False,
    #         "tickfont": dict(
    #             size=tick_lab_size,
    #             ),
    #         },
    #
    #     "xaxis": {
    #         "title": "Reaction Coordinate",
    #         "zeroline": True,
    #         "titlefont": dict(size=axes_lab_size),
    #         "showgrid": False,
    #
    #         # "showticklabels": False,
    #
    #         "ticktext": xax_labels,
    #         "tickvals": [1.5 * i + 0.5 for i in range(len(xax_labels))],
    #
    #         "tickfont": dict(
    #             size=tick_lab_size,
    #             ),
    #         },
    #     #__| -------------------------------------------------------------------
    #
    #     #| - Legend ------------------------------------------------------------
    #     "legend": {
    #         "traceorder": "normal",
    #         "font": dict(size=legend_size)
    #         },
    #     #__| -------------------------------------------------------------------
    #
    #     #| - Plot Size
    #     "width": 200 * 4.,
    #     "height": 200 * 3.,
    #     #__|
    #
    #     }
    #__|

    layout = plotly_fed_layout(plot_title=plot_title)

    return(dat_lst, layout)

    #__|

    #__|
