#!/usr/bin/env python

"""ORR energetics classes and methods.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import numpy as np
import pandas as pd

import plotly.plotly as py
import plotly.graph_objs as go

pd.options.mode.chained_assignment = None

from orr_reaction.orr_series import ORR_Free_E_Series
#__|

class ORR_Free_E_Plot:
    """ORR free energy diagram class.

    ACTUALLY THIS IS GOING TO BE A GENERAL ORR/OER CLASS NOW!!!!!!!!!!!!!!!!!!!



    Development Notes:
        # TODO Should we consider the case where the bulk energy is not 0, and
        we have to normalize all of the species energies by it?
    """

    #| - ORR_Free_E_Plot ******************************************************

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
        system_properties=None,
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
            properties=system_properties,
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
        font_family="Courier New, monospace",  # "Courier New, monospace"
        plot_width=680,
        plot_height=510,

        annotation_size=12,
        ):
        """

        Development notes:
            Move all plot parameters to this method, since the user will call
                this method to obtain the layout.
        """
        #| - plotly_fed_layout

        if self.rxn_type == "ORR":
            # xax_labels = ["$O_{2}$", "$*OOH$", "$*O$", "$*OH$", "$H_{2}O$"]
            # xax_labels = ["O2", "*OOH", "*O", "*OH", "H2O"]

            xax_labels = [
                "O<sub>2</sub>",
                "*OOH",
                "*O",
                "*OH",
                "H<sub>2</sub>O",
                ]

        elif self.rxn_type == "OER":
            # xax_labels = ["$H_{2}O$", "$*OH$", "$*O$", "$*OOH$", "$O_{2}$"]
            xax_labels = ["H2O", "*OH", "*O", "*OOH", "O2"]

        print(axes_lab_size)
        print(tick_lab_size)

        layout = {
            "title": plot_title,

            "font": {
                # "family": "Courier New, monospace",
                "family": font_family,
                "size": plot_title_size,
                "color": "black",
                },


            #| - Axes ---------------------------------------------------------
            "yaxis": {
                "title": "Free Energy (eV)",
                # "title": "$\\Delta G (ev)$",

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
            "legend": {
                "traceorder": "normal",
                "font": dict(size=legend_size)
                },

            "showlegend": self.show_legend,

            #__| --------------------------------------------------------------

            #| - Plot Size
            # "width": 200 * 4.,
            # "height": 200 * 3.,
            #__|

            # "paper_bgcolor": 'rgba(0,0,0,0)',
            "plot_bgcolor": 'rgba(0,0,0,0)',

            # "width": 9. * 37.795275591,
            # "height": 9 * 37.795275591,

            "width": plot_width,
            "height": plot_height,

            }

        if self.show_H_e_pairs_annotations:
            annotations = self.H_e_pairs_annotations(font_size=annotation_size)

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














    def create_scaling_relations_plot(self,
        y_ax_spec,
        x_ax_spec="oh",
        smart_format_dict=None,
        x_range=[0, 1.5]
        ):
        """Return plotly data and layout objects for scaling relations.

        Args:
            y_ax_spec:
            x_ax_spec:
        """
        #| - create_scaling_relations_plot

        #| - Internal Methods
        # TODO Should put these in a more accesible place

        def create_smart_format_dict(property_dict, smart_format_dict):
            """Create smart format dictionary.

            Args:
                property_dict:
                smart_format_dict:
            """
            format_dict = {}
            for key_i, value_i in property_dict.items():
                for format_i in smart_format_dict:
                    if list(format_i[0])[0] == key_i:
                        if list(format_i[0].values())[0] == value_i:
                            format_dict.update(format_i[1])

            return(format_dict)

        def create_series_name(series):
            """
            """
            #| - create_series_name
            name_i = ""
            for key, value in series_i.properties.items():
                if key == "coverage":
                    continue

                name_i += str(key) + ": " + str(value) + " | "

            return(name_i)
            #__|

        def ooh_oh_scaling(E_OH):
            return(E_OH + 3.2)

        def o_oh_scaling(E_OH):
            return(2 * E_OH)

        def create_trace_i(
            x_energy,
            y_energy,
            smart_format_i
            ):
            """
            """
            #| - create_trace_i
            trace_i = go.Scatter(
                x=x_energy,
                y=y_energy,
                text=name_i,
                name=name_i,
                mode='markers',
                marker=dict(
                    size=14,
                    color=smart_format_i["color2"],
                    line=dict(
                        color=smart_format_i["color1"],
                        width=4,
                        )
                    )
                )

            return(trace_i)
            #__|

        def create_layout(
            y_ax_spec,
            x_ax_spec,
            title="Scaling Relations",
            ):
            """
            """
            #| - create_layout
            if y_ax_spec == "ooh":
                y_ax_title = "G_OOH"
            elif y_ax_spec == "o":
                y_ax_title = "G_O"

            if x_ax_spec == "oh":
                x_ax_title = "G_OH"

            layout_i = dict(
                title=title,
                xaxis=dict(
                    title=x_ax_title,
                    zeroline=False,
                    ),
                yaxis=dict(
                    title=y_ax_title,
                    zeroline=False,
                    ),
                legend=dict(
                    x=0.,
                    y=1.8,
                    font=dict(
                        size=10,
                        ),
                    ),
                )

            return(layout_i)
            #__|

        #__|

        #| - Default Smart Format Dict
        if smart_format_dict is None:
            print("No smart format given!")
            smart_format_dict = [
                [{"bulk_system": "IrO3"}, {"color2": "green"}],
                [{"bulk_system": "IrO2"}, {"color2": "yellow"}],

                [{"coverage_type": "o_covered"}, {"symbol": "s"}],
                [{"coverage_type": "h_covered"}, {"symbol": "^"}],

                [{"facet": "110"}, {"color1": "red"}],
                [{"facet": "211"}, {"color1": "green"}],
                [{"facet": "100"}, {"color1": "black"}],
                ]
        #__|

        assert (x_ax_spec == "oh"), "Only *OH as the x-axis is allowed now"

        #| - Processing Data Points
        data_ooh_oh = []
        data_o_oh = []
        for series_i in self.series_list:
            e_oh = series_i.energy_states_dict["oh"]
            e_ooh = series_i.energy_states_dict["ooh"]
            e_o = series_i.energy_states_dict["o"]


            # print(series_i.properties)
            # print(smart_format_dict)
            # print("________________________")


            smart_format_i = create_smart_format_dict(
                series_i.properties,
                smart_format_dict,
                )

            print(smart_format_i)


            name_i = create_series_name(series_i)

            trace_i = create_trace_i(e_oh, e_ooh, smart_format_i)
            data_ooh_oh.append(trace_i)

            trace_i = create_trace_i(e_oh, e_o, smart_format_i)
            data_o_oh.append(trace_i)
        #__|

        #| - Ideal Scaling Lines
        scaling_trace = go.Scatter(
            x=[x_range[0], x_range[1]],
            y=[ooh_oh_scaling(x_range[0]), ooh_oh_scaling(x_range[1])],
            name='OOH_OH Scaling',
            mode='lines',
            line=dict(
                color="black",
                width=1,
                ),
            )
        data_ooh_oh.append(scaling_trace)

        scaling_trace = go.Scatter(
            x=[x_range[0], x_range[1]],
            y=[o_oh_scaling(x_range[0]), o_oh_scaling(x_range[1])],
            name='O_OH Scaling',
            mode='lines',
            line=dict(
                color="black",
                width=1,
                ),
            )
        data_o_oh.append(scaling_trace)
        #__|

        #| - Plot Layout Settings
        layout_ooh_oh = create_layout(
            y_ax_spec,
            x_ax_spec,
            title="OOH vs OH Scaling",
            )

        layout_o_oh = create_layout(
            y_ax_spec,
            x_ax_spec,
            title="O vs OH Scaling",
            )
        #__|

        if y_ax_spec == "ooh":
            return(data_ooh_oh, layout_ooh_oh)
        elif y_ax_spec == "o":
            return(data_o_oh, layout_o_oh)
        #__|

    #__| **********************************************************************
