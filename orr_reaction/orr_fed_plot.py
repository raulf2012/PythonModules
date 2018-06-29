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

class ORR_Free_E_Plot:
    """ORR free energy diagram class.

    Development Notes:
        # TODO Should we consider the case where the bulk energy is not 0, and
        we have to normalize all of the species energies by it?
    """

    #| - ORR_Free_E_Plot *******************************************************

    def __init__(self,
        free_energy_df=None,
        ORR_Free_E_series_list=None,
        # system_properties=None,
        state_title="adsorbate",
        free_e_title="ads_e",

        num_states=5,

        smart_format=None,

        bias=0.,

        # opt_name=None,
        # properties=None,

        color_list=None,

        # i_cnt=0,

        hover_text_col=None,

        # plot_mode="all",
        # smart_format=None,

        # Plotting ************************************************************
        show_H_e_pairs_annotations=True,
        show_legend=True,
        rxn_type="ORR",  # ORR and OER
        ):
        """
        Input variables to class instance.

        Args:
            free_energy_df:
                Pandas dataframe containing the adsorbates as rows
                Required columns, adsorbate, free energy
            system_properties:
            state_title:
            free_e_title:

            rxn_type:
                ORR or OER
        """
        #| - __init__

        #| - Setting Instance Attributes
        self.fe_df = free_energy_df
        # self.sys_props = system_properties
        self.state_title = state_title
        self.fe_title = free_e_title
        self.num_states = num_states

        # self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
        # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]  #COMBAK Is this being used anymore

        self.bias = bias
        self.color_list = color_list
        self.hover_text_col = hover_text_col
        self.smart_format = smart_format

        self.show_H_e_pairs_annotations = show_H_e_pairs_annotations
        self.show_legend = show_legend

        # ***********************************
        self.plot_states_sep = 0.3
        self.plot_states_width = 1.

        self.rxn_type = rxn_type
        #__|


        if self.rxn_type == "ORR":
            self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
            self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]

        elif self.rxn_type == "OER":
            self.rxn_mech_states = ["bulk", "oh", "o", "ooh", "bulk"]
            self.ideal_energy = [0, 1.23, 2.46, 3.69, 4.92]


        self.rxn_x_coord_array = self.create_rxn_coord_array(
            self.num_states,
            spacing=self.plot_states_sep,
            step_size=self.plot_states_width,
            )


        self.mid_state_x_array = self.create_mid_state_x_array()
        # x_array_data = self.rxn_x_coord_array


        if ORR_Free_E_series_list is None:
            self.series_list = []
        else:
            self.series_list = ORR_Free_E_series_list

        #| - __old__
        # if free_energy_df is not None:
        #     self.add_bulk_entry()
        #     self.fill_missing_data()
        #
        #     self.num_of_states = len(self.fe_df) + 1  # bulk, OOH, O, OH, bulk
        #     self.energy_lst = self.rxn_energy_lst()
        #     self.num_of_elec = range(self.num_of_states)[::-1]
        #     self.overpotential = self.calc_overpotential()[0]
        #     self.limiting_step = self.calc_overpotential()[1]
        #     # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]
        #     self.energy_lst_h2o2 = self.rxn_energy_lst_h2o2()
        #     self.overpotential_h2o2 = self.calc_overpotential_h2o2()
        #__|

        #__|

    def ideal_ORR_series(self,
        ):
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


        self.add_series(
            df_ideal,
            plot_mode="full_lines",  # ##########
            opt_name="Ideal ORR Catalyst",
            smart_format=False,

            # state_title=self.state_title,
            # free_e_title=self.fe_title,
            # bias=self.bias,
            #
            # # opt_name=None,  # #######
            # # properties=opt_name,
            # # color_list=self.color_list,
            # # i_cnt=0,  # ##########
            # hover_text_col=self.hover_text_col,
            #
            # # smart_format=self.smart_format,
            )

        #__|

    def add_series(self,
        fe_df,
        plot_mode="all",
        opt_name=None,
        smart_format=True,
        overpotential_type="ORR",
        ):
        """
        """
        #| - add_series
        if smart_format:
            smart_format_i = self.smart_format
        else:
            smart_format_i = None

        # free_energy_df=None,
        # # system_properties=None,
        # state_title="adsorbate",
        # free_e_title="ads_e",
        #
        # bias=0.,
        # rxn_x_coord_array=None,
        # opt_name=None,
        # properties=None,
        # color_list=None,
        # i_cnt=0,
        # hover_text_col=None,
        # plot_mode="all",
        # smart_format=None,
        # overpotential_type="ORR",
        # rxn_type="ORR",

        ORR_Series = ORR_Free_E_Series(
            free_energy_df=fe_df,
            # system_properties=None,
            state_title=self.state_title,
            free_e_title=self.fe_title,
            bias=self.bias,
            rxn_x_coord_array=self.rxn_x_coord_array,
            opt_name=opt_name,  # #######

            # properties=opt_name,
            color_list=self.color_list,
            i_cnt=0,  # ##########
            hover_text_col=self.hover_text_col,
            plot_mode=plot_mode,  # ##########
            smart_format=smart_format_i,

            # overpotential_type=self.rxn_type,
            rxn_type=self.rxn_type,
            )

        self.series_list.append(ORR_Series)
        #__|

    def plotly_data(self):
        """
        """
        #| - plotly_data
        master_data_list = []
        for series_i in self.series_list:
            master_data_list += series_i.series_plot

        return(master_data_list)
        #__|

    def plotly_fed_layout(self,
        plot_title="FED",
        plot_title_size=18,
        tick_lab_size=16,
        axes_lab_size=18,
        legend_size=18,
        # font_family="Computer Modern"  # "Courier New, monospace"
        font_family="Courier New, monospace"  # "Courier New, monospace"
        ):
        """

        Development notes:
            Move all plot parameters to this method, since the user will call
                this method to obtain the layout.
        """
        #| - plotly_fed_layout

        if self.rxn_type == "ORR":
            # xax_labels = ["$O_{2}$", "$*OOH$", "$*O$", "$*OH$", "$H_{2}O$"]
            xax_labels = ["O2", "*OOH", "*O", "*OH", "H2O"]

            # xax_labels = ["O<sub>2</sub>", "*OOH", "*O", "*OH", "H<sub>2</sub>O"]

        elif self.rxn_type == "OER":
            # xax_labels = ["$H_{2}O$", "$*OH$", "$*O$", "$*OOH$", "$O_{2}$"]
            xax_labels = ["H2O", "*OH", "*O", "*OOH", "O2"]

        layout = {
            "title": plot_title,

            "font": {
                # "family": "Courier New, monospace",
                "family": font_family,
                "size": plot_title_size,
                "color": "black",
                },

            #| - Axes --------------------------------------------------------------
            "yaxis": {
                "title": "Free Energy (eV)",
                # "title": "$\\Delta G (ev)$",

                "zeroline": True,
                "showline": True,
                "mirror": 'ticks',
                "showgrid": False,

                "titlefont": dict(size=axes_lab_size),

                "tickfont": dict(
                    size=tick_lab_size,
                    ),

                "autotick": False,
                "ticks": 'inside',
                "tick0": 0,
                "dtick": 0.5,
                "ticklen": 8,
                "tickwidth": 2,
                # "tickcolor": '#000'

                },

            "xaxis": {
                # "title": "Reaction Coordinate",

                "zeroline": True,
                "showline": True,
                "mirror": 'ticks',
                "showgrid": False,

                # "zeroline": True,
                # "showgrid": False,
                "titlefont": dict(size=axes_lab_size),

                "showticklabels": True,

                "ticktext": xax_labels,
                # "tickvals": [1.5 * i + 0.5 for i in range(len(xax_labels))],
                # "tickvals": [self.plot_states_width * i + 0.5
                #     for i in range(len(xax_labels))],
                "tickvals": self.mid_state_x_array,

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

            "showlegend": self.show_legend,

            }

        if self.show_H_e_pairs_annotations:
            annotations = self.H_e_pairs_annotations()

            if "annotations" in list(layout):
                layout["annotations"] += annotations
            else:
                layout["annotations"] = annotations

        return(layout)
        #__|


    def max_y_value_per_step(self):
        """
        """
        #| - max_y_value_per_step
        fe_matrix = []
        for series_i in self.series_list:
            fe_matrix.append(
                np.array(series_i.energy_lst),
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
        # ann_font_size = 18
        # states_sep = self.plot_states_sep
        # states_width = self.plot_states_width


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

    #__| **********************************************************************
