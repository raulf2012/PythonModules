#!/usr/bin/env python

"""ORR FED plotting class.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import numpy as np
import pandas as pd
import copy

from sklearn.linear_model import LinearRegression

import plotly.graph_objs as go

pd.options.mode.chained_assignment = None

from oxr_reaction.oxr_series import ORR_Free_E_Series
from oxr_reaction.adsorbate_scaling import lim_U_i
#__|


# ███████ ███████         ██████  ██       ██████  ████████
# ██      ██              ██   ██ ██      ██    ██    ██
# █████   █████           ██████  ██      ██    ██    ██
# ██      ██              ██      ██      ██    ██    ██
# ██      ███████ ███████ ██      ███████  ██████     ██

class Free_Energy_Plot():
    """
    Development Notes:
        Take the FED methods out of the ORR_Free_E_Plot class and into this one
    """

    #| - Free_Energy_Plot *****************************************************

    def __init__(self,
        ORR_Free_E_Plot,
        bias=0.,
        opt_name=None,
        properties=None,
        color_list=None,
        i_cnt=0,
        hover_text_col=None,
        plot_mode="all",
        smart_format=None,
        rxn_type="ORR",
        show_legend=True,
        show_H_e_pairs_annotations=True,
        ):
        """
        Input variables to class instance.

        Args:
            ORR_Free_E_Plot:
            mode:
                "all", "ooh_vs_oh", "o_vs_oh"
        """
        #| - __init__
        self.bias = bias
        self.opt_name = opt_name
        self.properties = properties
        self.color_list = color_list
        self.i_cnt = i_cnt
        self.hover_text_col = hover_text_col
        self.plot_mode = plot_mode
        self.smart_format = smart_format
        self.rxn_type = rxn_type
        self.show_H_e_pairs_annotations = show_H_e_pairs_annotations

        self.ORR_Free_E_Plot = ORR_Free_E_Plot

        # #####################################################################
        self.plot_states_sep = 0.3
        self.plot_states_width = 1.

        self.show_legend = show_legend

        self.rxn_x_coord_array = self.create_rxn_coord_array(
            self.ORR_Free_E_Plot.num_states,
            spacing=self.plot_states_sep,
            step_size=self.plot_states_width,
            )

        self.mid_state_x_array = self.create_mid_state_x_array()
        #__|

    def create_fed_plot(self,

        ):
        """
        """
        #| - create_fed_plot
        plot_data = []
        for series_i in self.ORR_Free_E_Plot.series_list:

            data_i = self.plot_fed_series(
                series_i,

                bias=self.bias,
                opt_name=self.opt_name,
                properties=self.properties,
                color_list=self.color_list,
                i_cnt=self.i_cnt,
                hover_text_col=self.hover_text_col,
                plot_mode=self.plot_mode,
                smart_format=self.smart_format,
                overpotential_type=self.rxn_type,
                )

            plot_data += data_i

        return(plot_data)
        #__|

    def ideal_ORR_series(self):
        """
        """
        #| - ideal_ORR_series
        # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]

        ideal_data_list = [

            {
                "adsorbate": "ooh",
                "ads_e": 3.69,
                },

            {
                "adsorbate": "o",
                "ads_e": 2.46,
                },

            {
                "adsorbate": "oh",
                "ads_e": 1.23,
                },

            ]

        df_ideal = pd.DataFrame(ideal_data_list)

        # system_properties={"system": sys_name_i},
        # name_i=sys_name_i,
        # format_dict=format_dict_i,
        # opt_name=sys_name_i.replace("_", " "),

        self.ORR_Free_E_Plot.add_series(
            df_ideal,
            plot_mode="full_lines",  # ##########
            format_dict={"opacity": 0.},
            name_i="Ideal ORR Catalyst",
            opt_name="Ideal ORR Catalyst",
            smart_format=False,
            color=None,
            )

        #__|

    def create_rxn_coord_array(self,
        rxn_steps,
        spacing=0,
        step_size=1,
        ):
        """
        Create a reaction coordinate array ([0, 1, 1, 2, 2, 3]) for plotting.

        Args:
            rxn_steps: <type 'int'>
                Number of steps in reaction coordinate including initial
                and final state.
                Ex. A -> Intermediate -> C has 3 steps/states

            spacing: <type 'float'>
                Spacing inbetween the energy levels. The default of 0 creates
                a free energy diagram that looks like steps
        """
        #| - create_rxn_coord_array
        lst = []
        for i in range(1, rxn_steps):
            if i == 1:
                lst.append(step_size)
                lst.append(step_size + spacing)
            if i != 1:
                lst.append(lst[-1] + step_size)
                lst.append(lst[-2] + step_size + spacing)

        lst.insert(0, 0)
        lst.append(lst[-1] + step_size)

        return(lst)
        #__|


    def plot_fed_series(self,
        series_i,

        bias=0.,
        opt_name=None,
        properties=None,
        color_list=None,
        i_cnt=0,
        hover_text_col=None,
        plot_mode="all",
        smart_format=None,
        overpotential_type="ORR",
        ):
        """
        Process data for FED plot.

        Args:
            bias:
            opt_name
            properties:
            color_list:
            i_cnt:
            hover_text_col:
                Dataframe column name to be used for hover text

        #FIXME | This is  fairly rough as of right now
        """
        #| - plot_fed_series
        e_list = series_i.energy_lst
        e_list = series_i.apply_bias(bias, e_list)

        for n, i in enumerate(e_list):
            if np.isnan(i) is True:
                e_list[n] = None

        if color_list is None:
            color_list = ["red"]

        # name_i = self.series_name
        name_i = series_i.series_name

        #| - Hover Text
        if hover_text_col is not None:
            if type(hover_text_col) is not list:
                hover_text_list = self.property_list(hover_text_col)

            else:
                hover_text_lists = []
                for col_i in hover_text_col:

                    # replacing nan with ""
                    tmp = self.property_list(col_i)
                    hover_text_i = ["" if x is np.nan else x for x in tmp]
                    hover_text_lists.append(hover_text_i)

                # TODO Removed last trailing " | "
                hover_text_list = []
                for items in zip(*hover_text_lists):

                    if all([True if i == "" else False for i in items]) is True:
                        hover_col_state_i = ""
                    else:
                        hover_col_state_i = " | ".join(items)

                    hover_text_list.append(hover_col_state_i)

        else:
            hover_text_list = [np.nan for j_cnt in list(range(5))]
        #__|

        #| - TEMP Picking color from "color_list" or "color" variable
        if series_i.color is not None:
            color_i = self.color
        else:
            color_i = color_list[i_cnt - 1]
        #__|

        dat_lst = self.__create_plotly_series__(
            series_i,
            e_list,
            group=name_i,
            name=name_i,
            hover_text=hover_text_list,
            # color=color_list[i_cnt - 1],
            color=color_i,
            plot_mode=plot_mode,
            smart_format=smart_format,
            )

        return(dat_lst)

        #__|


    def __create_plotly_series__(self,
        series_i,
        energy_lst,
        name="TEMP",
        group="group1",
        hover_text=None,
        color="rgb(22, 96, 167)",
        plot_mode="all",
        smart_format=None,
        ):
        """
        Create a plotly series for the current instance.

        Args:
            energy_lst:
            name:
            group:
            color:
            plot_mode:
                "all"
                "states_only"
                "full_lines"
        """
        #| - create_plotly_series
        y_dat = self.__convert_to_plotting_list__(energy_lst)

        if hover_text is None:
            hover_text = [np.nan for i_ind in range(5)]


        #| - Parameters
        if plot_mode == "all":
            show_leg_2 = False
        elif plot_mode == "states_only":
            show_leg_2 = False
        elif plot_mode == "full_lines":
            show_leg_2 = True
        #__|

        #| - Adding Breaks in Data
        x_dat = self.rxn_x_coord_array

        new_x_dat = copy.copy(x_dat)
        new_y_dat = copy.copy(y_dat)

        cnt = 2
        for i_ind in range(int(len(x_dat) / 2 - 1)):
            fill = new_x_dat[cnt - 1]
            new_x_dat.insert(cnt, fill)
            new_y_dat.insert(cnt, None)
            cnt += 3
        #__|

        #| - Creating x-data in middle of states
        short_y = np.array(y_dat)[::2]

        xdat = list(set(new_x_dat))
        xdat.sort()

        cnt = 0
        short_x = []
        for i_ind in range(int(len(xdat) / 2)):
            short_x.append(xdat[cnt] + 0.5)  # FIXME Replace 0.5 with variable
            cnt += 2
        #__|

        #| - Smart Format Dict ************************************************

        #| - DICTS
        plot_parameter_dict = {
            "dash": None,
            }

        # smart_format = [
        #
        #     # [
        #     #     {"spinpol": True},
        #     #     {"dash": "dot"},
        #     #     ],
        #
        #     [
        #         {"system": "Fe_slab"},
        #         {"dash": "dashdot"},
        #         ],
        #
        #     [
        #         {"system": "N_graph_Fe"},
        #         {"dash": "dot"},
        #         ],
        #
        #     [
        #         {"system": "graph_Fe"},
        #         {"dash": "dash"},
        #         ],
        #
        #     [
        #         {"system": "graphene"},
        #         {"dash": None},
        #         ],
        #
        #     ]

        #__|

        # if self.fe_df is not None and smart_format is not None:
        if series_i.fe_df is not None and smart_format is not None:
            for format_i in smart_format:

                # format_i key:values
                df_col_name = list(format_i[0])[0]
                value_i = list(format_i[0].values())[0]
                setting_name = list(format_i[1])[0]
                setting_value = list(format_i[1].values())[0]


                if df_col_name in list(self.fe_df):

                    # *OOH, *O, *OH entries must have value_i as the
                    # column entries in the column df_col_name

                    df = self.fe_df
                    df_wo_bulk = df[df["adsorbate"] != "bulk"]

                    if all(df_wo_bulk[df_col_name] == value_i):
                        plot_parameter_dict.update(
                            {setting_name: setting_value},
                            )

                else:
                    print("Dataframe column " + df_col_name + " not present")

        #__| ******************************************************************

        #| - Series Color
        if "color" in list(plot_parameter_dict):
            color_out = plot_parameter_dict["color"]
        else:
            plot_parameter_dict["color"] = color
        #__|

        if series_i.format_dict:
            format_i = series_i.format_dict
        else:
            format_i = plot_parameter_dict

        #| - Plotly Scatter Plot Instances ************************************

        #| - Thick horizontal state lines
        data_1 = go.Scatter(
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            showlegend=False,
            name=name,
            hoverinfo="none",  # TEMP - 180317

            # text=hover_text,

            connectgaps=False,
            line=dict(
                # color=color,
                # color=plot_parameter_dict["color"],
                # color=format_i["color_2"],
                color=format_i.get("color_2", "red"),

                # COMBAK
                # width=2,
                width=4,

                # dash="dot",  # TEMP
                dash=plot_parameter_dict["dash"],  # TEMP
                ),
            mode="lines",
            )
        #__|

        #| - Full, thin line
        data_2 = go.Scatter(
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            name=name,
            connectgaps=True,
            showlegend=show_leg_2,
            hoverinfo="none",
            text=hover_text,
            line=dict(
                color=format_i.get("color_2", "red"),
                width=1,
                ),
            mode="lines",
            )
        #__|

        #| - Points in middle of energy states (For convienient hover)

        data_3 = go.Scatter(
            x=short_x,
            y=short_y,
            legendgroup=group,
            name=name,
            showlegend=True,
            hoverinfo="y+text",
            # hoverinfo="y+name",
            # text=hover_text,
            text=name,
            marker=dict(
                size=format_i.get("marker_size", 0),
                color=format_i.get("color_2", "black"),
                # opacity=0.,
                # opacity=0.8,
                opacity=format_i.get("opacity", 0.8),

                symbol=format_i.get("symbol_i", "circle"),
                line=dict(
                    # color=smart_format_i[marker_border_color_key],
                    color=format_i.get("color_1", "black"),
                    width=1.,
                    )

                ),
            mode="markers",
            )
        #__|

        #| - Points in middle of RDS states
        # HACK
        if series_i.limiting_step == ["ooh", "bulk"]:
            ind_i = 3
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]

        if series_i.limiting_step == ["bulk", "ooh"]:
            ind_i = 0
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]



        if series_i.limiting_step == ["o", "ooh"]:
            ind_i = 2
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]

        if series_i.limiting_step == ["ooh", "o"]:
            ind_i = 1
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]



        if series_i.limiting_step == ["oh", "o"]:
            ind_i = 1
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]

        if series_i.limiting_step == ["o", "oh"]:
            ind_i = 2
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]



        if series_i.limiting_step == ["bulk", "oh"]:
            ind_i = 0
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]

        if series_i.limiting_step == ["oh", "bulk"]:
            ind_i = 3
            short_x = short_x[ind_i:ind_i + 2]
            short_y = short_y[ind_i:ind_i + 2]

        # TEMP
        # print(series_i.limiting_step)
        # print(short_x)
        # print(short_y)
        # print("***********************************************")

        rds_data = go.Scatter(
            x=short_x,
            y=short_y,
            legendgroup=group,
            name=name,
            showlegend=False,
            hoverinfo="none",
            # text=hover_text,
            marker=dict(
                size=format_i.get("marker_size", 10) + 14,
                # color=format_i.get("color_2", "black"),
                color="grey",
                opacity=0.3,
                # symbol=format_i.get("symbol_i", "circle"),
                symbol="circle",
                line=dict(
                    # color=smart_format_i[marker_border_color_key],
                    color=format_i.get("color_1", "black"),
                    width=0.,
                    )

                ),
            mode="markers",
            )
        #__|

        #__| ******************************************************************

        #| - Plot Mode (which data series to plot)
        if plot_mode == "all":
            data_lst = [
                rds_data,
                data_1,
                data_2,
                data_3,
                ]
        elif plot_mode == "states_only":
            data_lst = [
                rds_data,
                data_1,
                data_3,
                ]
        elif plot_mode == "full_lines":
            data_lst = [
                rds_data,
                data_2,
                data_3,
                ]
        #__|

        return(data_lst)
        #__|

    def __convert_to_plotting_list__(self,
        energy_lst,
        spacing=0.5,
        step_size=1,
        ):
        """Repeat entries in energy list to conform to FED plot.

        Modifies an energy list for plotting by repeating each entry
        Ex. [4.92, 3.69, ... ] -> [4.92, 4.92, 3.69, 3.69, ... ]

        Args:
            energy_lst: <type 'list'>
            spacing:
            step_size:
        """
        #| - __convert_to_plotting_list__
        tmp_list = range(len(energy_lst) * 2)
        energy_dupl_lst = [energy_lst[i // 2] for i in tmp_list]

        # rxn_coord_steps = self.create_rxn_coord_array(
        #     len(energy_lst),
        #     spacing=spacing,
        #     step_size=step_size,
        #     )
        # out_list = [rxn_coord_steps, energy_dupl_lst]

        return(energy_dupl_lst)
        #__|

    def max_y_value_per_step(self):
        """
        """
        #| - max_y_value_per_step
        fe_matrix = []

        # for series_i in self.series_list:
        for series_i in self.ORR_Free_E_Plot.series_list:

            energy_list_i = series_i.apply_bias(self.bias, series_i.energy_lst)

            fe_matrix.append(
                energy_list_i,
                # np.array(series_i.energy_lst),
                )
        fe_matrix = np.array(fe_matrix)

        max_y_val_list = []
        for step_i in range(fe_matrix.shape[1]):
            max_y_val_list.append(fe_matrix[:, step_i].max())

        return(max_y_val_list)
        #__|


    def H_e_pairs_annotations(self,
        font_size=18,
        ):
        """

        Args:
            font_size:
        """
        #| - H_e_pairs_annotations
        mid_state_x_array = self.mid_state_x_array

        rxn_x_array = self.rxn_x_coord_array
        max_y_val_list = self.max_y_value_per_step()

        def add_annot(
            ind,
            rxn_x_array,
            step,
            annotations,
            y=5.5,
            font_size=18,
            text="TEMP",
            mid_state_x_array=None,
            ):
            """Append annotation object to annotations list.

            Args:
                ind:
                rxn_x_array:
                step:
                annotations:
                y:
                font_size:
                text:
            """
            #| - add_annot
            ann_i = dict(
                # x=(rxn_x_array[ind] + rxn_x_array[ind + step]) / 2.,
                x=mid_state_x_array[ind],
                y=y + 0.4,
                xref='x',
                yref='y',
                text=text,
                showarrow=False,
                font=dict(
                    color="black",
                    size=font_size
                    ),
                ),
            annotations += ann_i
            #__|

        annotations = []

        add_annot(
            0, rxn_x_array, 1, annotations,
            y=max_y_val_list[0], font_size=font_size,
            text="4(H<sup>+</sup> + e<sup>-</sup>)",
            mid_state_x_array=mid_state_x_array,
            )

        add_annot(
            1, rxn_x_array, 2, annotations,
            y=max_y_val_list[1], font_size=font_size,
            text="3(H<sup>+</sup> + e<sup>-</sup>)",
            mid_state_x_array=mid_state_x_array,
            )

        add_annot(
            # 3, rxn_x_array, 2, annotations,
            2, rxn_x_array, 2, annotations,
            y=max_y_val_list[2], font_size=font_size,
            text="2(H<sup>+</sup> + e<sup>-</sup>)",
            mid_state_x_array=mid_state_x_array,
            )

        add_annot(
            # 5, rxn_x_array, 2, annotations,
            3, rxn_x_array, 2, annotations,
            y=max_y_val_list[3], font_size=font_size,
            text="1(H<sup>+</sup> + e<sup>-</sup>)",
            mid_state_x_array=mid_state_x_array,
            )

        add_annot(
            # 7, rxn_x_array, 2, annotations,
            4, rxn_x_array, 2, annotations,
            y=max_y_val_list[4], font_size=font_size,
            text="0(H<sup>+</sup> + e<sup>-</sup>)",
            mid_state_x_array=mid_state_x_array,
            )

        return(annotations)
        #__|


    def create_mid_state_x_array(self):
        """
        """
        #| - create_mid_state_x_array
        x_array_data = self.rxn_x_coord_array
        state_width = self.plot_states_width

        xdat = list(set(x_array_data))
        xdat.sort()

        cnt = 0
        short_x = []
        for i_ind in range(int(len(xdat) / 2)):
            short_x.append(xdat[cnt] + state_width / 2.)
            cnt += 2

        return(short_x)
        #__|


    def get_plotly_layout(self, layout_dict=None):
        """
        """
        #| - plotly_fed_layout
        tick_lab_size = 16
        axes_lab_size = 18
        legend_size=18
        annotation_size = 12

        #| - OER vs ORR Settings
        if self.rxn_type == "ORR":
            xax_labels = [
                "O<sub>2</sub>",
                "*OOH",
                "*O",
                "*OH",
                "H<sub>2</sub>O",
                ]

        elif self.rxn_type == "OER":
            # xax_labels = ["$H_{2}O$", "$*OH$", "$*O$", "$*OOH$", "$O_{2}$"]
            xax_labels = [
                "H<sub>2</sub>O",
                "*OH",
                "*O",
                "*OOH",
                "O<sub>2</sub>",
                ]
        #__|

        #| - Layout
        layout = {

            "title": None,

            #| - Font Settings
            "font": {
                "family": "Arial",  # "Courier New, monospace"
                "color": "black",
                },
            #__|

            #| - Margins ------------------------------------------------------
            "margin": go.layout.Margin(
                b=50.,
                l=50.,
                r=30.,
                t=30.,
                ),
            #__|

            #| - Axes ---------------------------------------------------------
            "yaxis": {
                "title": "Free Energy (eV)",
                "zeroline": False,
                "linecolor": 'black',
                "showline": True,
                "mirror": 'ticks',
                "showgrid": False,
                "titlefont": dict(size=axes_lab_size),
                "tickfont": dict(
                    size=tick_lab_size,
                    ),

                # "autotick": False,
                "ticks": 'inside',
                "tick0": 0,
                "dtick": 1.0,
                "ticklen": 2,
                "tickwidth": 1,
                "tickcolor": 'black',
                },

            "xaxis": {
                "title": "Reaction Coordinate",

                "zeroline": False,
                "linecolor": 'black',
                "showline": True,
                "mirror": 'ticks',
                "showgrid": False,

                # "zeroline": True,
                # "showgrid": False,
                "titlefont": dict(size=axes_lab_size),

                "showticklabels": True,

                "ticks": "",
                "ticktext": xax_labels,
                # "tickvals": [1.5 * i + 0.5 for i in range(len(xax_labels))],
                # "tickvals": [self.plot_states_width * i + 0.5
                #     for i in range(len(xax_labels))],
                "tickvals": self.mid_state_x_array,
                "tickcolor": 'black',
                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                },
            #__| --------------------------------------------------------------

            #| - Legend -------------------------------------------------------

            "legend": go.layout.Legend(
                x=1.1,
                xanchor=None,
                y=1.,
                yanchor="top",
                font=dict(size=legend_size),
                bgcolor="rgba(0,0,0,0.01)",

                arg=None,
                bordercolor=None,
                borderwidth=None,
                itemclick=None,
                itemdoubleclick=None,
                itemsizing=None,
                orientation=None,
                tracegroupgap=None,
                traceorder=None,
                uirevision=None,
                valign=None,
                ),

            "showlegend": self.show_legend,

            #__| --------------------------------------------------------------

            "paper_bgcolor": 'rgba(240,240,240,0.9)',

            }

        #__|

        #| - H/e Count Annotations
        if self.show_H_e_pairs_annotations:
            annotations = self.H_e_pairs_annotations(font_size=annotation_size)
            if "annotations" in list(layout):
                layout["annotations"] += annotations
            else:
                layout["annotations"] = annotations
        #__|

        layout = go.Layout(**layout)

        if layout_dict is not None:
            layout.update(layout_dict)

        return(layout)
        #__|

    #__| **********************************************************************
