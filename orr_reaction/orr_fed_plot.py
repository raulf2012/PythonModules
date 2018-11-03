#!/usr/bin/env python

"""ORR energetics classes and methods.

Author: Raul A. Flores
"""


#| - IMPORT MODULES
import numpy as np
import pandas as pd

# import itertools

from sklearn.linear_model import LinearRegression

# import plotly.plotly as py
import plotly.graph_objs as go

pd.options.mode.chained_assignment = None

from orr_reaction.orr_series import ORR_Free_E_Series
from orr_reaction.adsorbate_scaling import lim_U_i
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
            color=None,
            )

        #| - __old__

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

        #__|

        #__|

    def __create_series_name__(self, series_i):
        """
        """
        #| - __create_series_name__

        if series_i.properties is not None:
            name_i = ""
            for key, value in series_i.properties.items():

                # TODO | This is an old feature, get rid of it
                if key == "coverage":
                    continue

                name_i += str(key) + ": " + str(value) + " | "

        else:
            name_i = ""

        return(name_i)
        #__|

    def __create_smart_format_dict__(self, property_dict, smart_format_dict):
        """Create smart format dictionary.

        Args:
            property_dict:
            smart_format_dict:
        """
        #| - __create_smart_format_dict__
        if property_dict is None:
            return({})

        format_dict = {}
        for key_i, value_i in property_dict.items():
            for format_i in smart_format_dict:
                if list(format_i[0])[0] == key_i:
                    if list(format_i[0].values())[0] == value_i:
                        format_dict.update(format_i[1])

        return(format_dict)
        #__|

    def add_series(self,
        fe_df,
        plot_mode="all",
        name_i=None,
        group=None,
        opt_name=None,
        smart_format=True,
        overpotential_type="ORR",
        system_properties=None,
        property_key_list=None,
        color=None,
        ):
        """Add ORR_Free_E_Series instance to ORR_Free_E_Plot.series_list.

        Note: It would be much better to simply take all of the
        ORR_Free_E_Series arguments as a **kwargs term.

        Args:
            TEMP
        """
        #| - add_series
        if smart_format:
            smart_format_i = self.smart_format
        else:
            smart_format_i = None

        #| - __old__
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

        # print(system_properties)
        #__|

        ORR_Series = ORR_Free_E_Series(
            free_energy_df=fe_df,
            properties=system_properties,
            property_key_list=property_key_list,
            state_title=self.state_title,
            free_e_title=self.fe_title,
            group=group,
            bias=self.bias,
            rxn_x_coord_array=self.rxn_x_coord_array,
            name_i=name_i,
            opt_name=opt_name,  # #######

            # properties=opt_name,
            color_list=self.color_list,
            color=color,
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
                "font": dict(size=legend_size),
                "x": -0.1,
                "y": -1.2,
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


    # Deprecated **************************************************************
    def create_scaling_relations_plot(self,
        y_ax_spec,
        x_ax_spec="oh",
        smart_format_dict=None,

        x_range_ooh_vs_oh=[-1., 3.5],
        y_range_ooh_vs_oh=[0., 5.],
        x_range_o_vs_oh=[-1., 3.5],
        y_range_o_vs_oh=[0., 5.],
        x_range_oh_vs_oh=[-1., 4.],
        y_range_oh_vs_oh=[-1., 4.],
        ):
        """Return plotly data and layout objects for scaling relations.

        Args:
            y_ax_spec:
            x_ax_spec:
        """
        #| - create_scaling_relations_plot

        print("#########################################")
        print("DEPRECATED!!!!!!!!!!!!!")
        print("Use the new class Scaling_Relations_Plot")
        print("#########################################")

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
            """Return the *OOH adsorption energy given DG_*OH by scaling.

            Args:
                E_OH:DG_*OH energy of adsorption
            """
            #| - ooh_oh_scaling
            return(E_OH + 3.2)

            #__|

        def o_oh_scaling(E_OH):
            """Return the *OOH adsorption energy given DG_*OH by scaling.

            Args:
                E_OH: DG_*OH energy of adsorption.
            """
            #| - o_oh_scaling
            return(2 * E_OH)
            #__|

        def oh_oh_scaling(E_OH):
            """Return the *OH adsorption energy given DG_*OH by scaling.

            NOTE: TRIVIAL QUANTITY!!!!!!!!!!!!!!!!!!!

            Args:
                E_OH: DG_*OH energy of adsorption.
            """
            #| - oh_oh_scaling
            return(E_OH)
            #__|


        def create_trace_i(
            x_energy,
            y_energy,
            smart_format_i
            ):
            """
            """
            #| - create_trace_i
            # NOTE Looks like I need to put these in a list here
            x_energy = [x_energy]
            y_energy = [y_energy]

            trace_i = go.Scatter(
                x=x_energy,
                y=y_energy,
                text=name_i,
                name=name_i,
                mode='markers',
                marker=dict(
                    size=14,
                    symbol=smart_format_i.get("symbol", "circle"),
                    color=smart_format_i.get("color2", "pink"),
                    line=dict(
                        # color=smart_format_i["color1"],
                        color=smart_format_i.get("color1", "black"),
                        width=2,
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
                y_ax_title = "G<sub>ads,*OOH</sub> (eV)"
            elif y_ax_spec == "o":
                y_ax_title = "G<sub>ads,*O</sub> (eV)"
            elif y_ax_spec == "oh":
                y_ax_title = "G<sub>ads,*OH</sub> (eV)"


            if x_ax_spec == "oh":
                x_ax_title = "G<sub>ads,*OH</sub> (eV)"
            else:
                print("Only 'oh' is supported as the x-axis variable")

            tick_lab_size = 12 * (4. / 3.)
            axes_lab_size = 14 * (4. / 3.)
            # legend_size = 18

            #| - Common Axis Dict
            common_axis_dict = {

                # "range": y_axis_range,
                "zeroline": False,
                "showline": True,
                "mirror": 'ticks',
                "linecolor": 'black',
                "showgrid": False,

                "titlefont": dict(size=axes_lab_size),
                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                "ticks": 'inside',
                "tick0": 0,
                "tickcolor": 'black',
                # "dtick": 0.25,
                "ticklen": 2,
                "tickwidth": 1,
                }
            #__|

            # x_range_ooh_vs_oh=[0., 3.5],
            # y_range_ooh_vs_oh=[0., 5.],
            # x_range_o_vs_oh=[0., 3.5],
            # y_range_o_vs_oh=[0., 5.],

            if y_ax_spec == "ooh":
                x_range = x_range_ooh_vs_oh
            elif y_ax_spec == "o":
                x_range = x_range_o_vs_oh
            elif y_ax_spec == "oh":
                x_range = x_range_oh_vs_oh
            else:
                print("Woops - create_layout")

            if y_ax_spec == "ooh":
                y_range = y_range_ooh_vs_oh
            elif y_ax_spec == "o":
                y_range = y_range_o_vs_oh
            elif y_ax_spec == "oh":
                y_range = y_range_oh_vs_oh
            else:
                print("Woops - create_layout")


            layout_i = {
                "title": title,
                "titlefont": go.layout.Titlefont(size=24),

                "xaxis": dict(
                    common_axis_dict,
                    **{
                        "title": x_ax_title,
                        "range": x_range,
                        },
                    ),

                "yaxis": dict(
                    common_axis_dict,
                    **{
                        "title": y_ax_title,
                        "range": y_range,
                        },
                    ),

                "font": dict(
                    family='Arial',
                    # size=18,
                    color='black',
                    ),

                "width": 1.5 * 18.7 * 37.795275591,
                "height": 18.7 * 37.795275591,

                "showlegend": True,

                "legend": dict(
                    # x=0.,
                    # y=1.8,
                    font=dict(
                        size=10,
                        ),
                    ),
                }

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
        data_oh_oh = []
        for series_i in self.series_list:
            e_oh = series_i.energy_states_dict["oh"]
            e_ooh = series_i.energy_states_dict["ooh"]
            e_o = series_i.energy_states_dict["o"]

            smart_format_i = create_smart_format_dict(
                series_i.properties,
                smart_format_dict,
                )

            name_i = create_series_name(series_i)

            if series_i.color is not None:
                smart_format_i["color2"] = series_i.color

            trace_i = create_trace_i(e_oh, e_ooh, smart_format_i)
            data_ooh_oh.append(trace_i)

            trace_i = create_trace_i(e_oh, e_o, smart_format_i)
            data_o_oh.append(trace_i)

            trace_i = create_trace_i(e_oh, e_oh, smart_format_i)
            data_oh_oh.append(trace_i)

        #__|

        #| - Ideal Scaling Lines
        scaling_trace = go.Scatter(
            x=[x_range_ooh_vs_oh[0], x_range_ooh_vs_oh[1]],
            y=[
                ooh_oh_scaling(x_range_ooh_vs_oh[0]),
                ooh_oh_scaling(x_range_ooh_vs_oh[1]),
                ],
            name='*OOH vs *OH Scaling',
            mode='lines',
            line=dict(
                color="black",
                width=1,
                ),
            )
        data_ooh_oh.append(scaling_trace)

        scaling_trace = go.Scatter(
            x=[x_range_o_vs_oh[0], x_range_o_vs_oh[1]],
            y=[
                o_oh_scaling(x_range_o_vs_oh[0]),
                o_oh_scaling(x_range_o_vs_oh[1]),
                ],
            name='*O vs *OH Scaling',
            mode='lines',
            line=dict(
                color="black",
                width=1,
                ),
            )
        data_o_oh.append(scaling_trace)

        scaling_trace = go.Scatter(
            x=[x_range_oh_vs_oh[0], x_range_oh_vs_oh[1]],
            y=[
                oh_oh_scaling(x_range_oh_vs_oh[0]),
                oh_oh_scaling(x_range_oh_vs_oh[1]),
                ],
            name='*OH vs *OH Scaling',
            mode='lines',
            line=dict(
                color="black",
                width=1,
                ),
            )
        data_oh_oh.append(scaling_trace)

        #__|

        #| - Plot Layout Settings

        layout_ooh_oh = create_layout(
            y_ax_spec,
            x_ax_spec,
            title="*OOH vs *OH Scaling",
            )

        layout_o_oh = create_layout(
            y_ax_spec,
            x_ax_spec,
            title="*O vs *OH Scaling",
            )

        layout_oh_oh = create_layout(
            y_ax_spec,
            x_ax_spec,
            title="*OH vs *OH Scaling",
            )
        #__|

        if y_ax_spec == "ooh":
            return(data_ooh_oh, layout_ooh_oh)
        elif y_ax_spec == "o":
            return(data_o_oh, layout_o_oh)
        elif y_ax_spec == "oh":
            return(data_oh_oh, layout_oh_oh)
        #__|

    #__| **********************************************************************


class Scaling_Relations_Plot():
    """Plot scaling relations and some simple fitting schemes.

    Development Notes:
        TEMP
    """

    #| - Scaling_Relations_Plot ***********************************************

    def __init__(self,
        ORR_Free_E_Plot,
        mode="all",

        plot_range={
            "y": [1., 5.],
            "x": [-2., 4.],
            },

        x_ax_species="oh",
        ):
        """
        Input variables to class instance.

        Args:
            ORR_Free_E_Plot:
            mode:
                "all", "ooh_vs_oh", "o_vs_oh"
        """
        #| - __init__
        self.ORR_Free_E_Plot = ORR_Free_E_Plot

        assert (x_ax_species == "oh"), "Only *OH as the x-axis is allowed now"


        self.data_points = {
            "ooh_vs_oh": [],
            "o_vs_oh": [],
            "oh_vs_oh": [],
            }
        self.data_lines = []

        self.x_range = plot_range["x"]
        self.y_range = plot_range["y"]

        # self.layout = self.__create_layout__(
        #     title="Scaling Relations",
        #     showlegend=True,
        #     )


        #__|

    def create_scaling_relations_plot(self,
        smart_format_dict=None,
        ):
        """Return plotly data and layout objects for scaling relations.

        Args:
            y_ax_spec:
            x_ax_spec:
        """
        #| - create_scaling_relations_plot

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

        #| - Processing Data Points
        for series_i in self.ORR_Free_E_Plot.series_list:

            e_oh = series_i.energy_states_dict["oh"]
            e_ooh = series_i.energy_states_dict["ooh"]
            e_o = series_i.energy_states_dict["o"]


            # Change self.__create_smart_format_dict__ to,
            # self.ORR_Free_E_Plot.__create_smart_format_dict__
            # TODO

            smart_format_i = self.ORR_Free_E_Plot.__create_smart_format_dict__(
                series_i.properties,
                smart_format_dict,
                )

            # This is the old way of getting the name
            # name_i = self.ORR_Free_E_Plot.__create_series_name__(series_i)

            name_i = series_i.series_name

            if series_i.color is not None:
                smart_format_i["color2"] = series_i.color

            #| - ooh_vs_oh
            trace_i = self.__create_trace_i__(
                e_oh,
                e_ooh,
                smart_format_i,
                name_i,
                legendgroup="ooh_vs_oh",
                )
            # self.data_ooh_oh.append(trace_i)
            self.data_points["ooh_vs_oh"].append(trace_i)
            #__|

            #| - o_vs_oh
            trace_i = self.__create_trace_i__(
                e_oh,
                e_o,
                smart_format_i,
                name_i,
                legendgroup="o_vs_oh",
                )
            # self.data_o_oh.append(trace_i)
            self.data_points["o_vs_oh"].append(trace_i)
            #__|

            #| - oh_vs_oh
            trace_i = self.__create_trace_i__(
                e_oh,
                e_oh,
                smart_format_i,
                name_i,
                legendgroup="oh_vs_oh",
                )
            # self.data_oh_oh.append(trace_i)
            self.data_points["oh_vs_oh"].append(trace_i)
            #__|

        #__|

        #__|

    # Deprecated, delete this later
    def __create_smart_format_dict__(self, property_dict, smart_format_dict):
        """Create smart format dictionary.

        Args:
            property_dict:
            smart_format_dict:
        """
        #| - __create_smart_format_dict__
        format_dict = {}
        for key_i, value_i in property_dict.items():
            for format_i in smart_format_dict:
                if list(format_i[0])[0] == key_i:
                    if list(format_i[0].values())[0] == value_i:
                        format_dict.update(format_i[1])

        return(format_dict)
        #__|

    def __create_series_name__(self, series_i):
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

    def __create_trace_i__(self,
        x_energy,
        y_energy,
        smart_format_i,
        name_i,
        legendgroup=None,
        ):
        """
        """
        #| - create_trace_i
        # NOTE Looks like I need to put these in a list here
        x_energy = [x_energy]
        y_energy = [y_energy]

        trace_i = go.Scatter(
            x=x_energy,
            y=y_energy,
            text=name_i,
            name=name_i,
            mode="markers",
            legendgroup=legendgroup,
            marker=dict(
                size=14,
                symbol=smart_format_i.get("symbol", "circle"),
                color=smart_format_i.get("color2", "pink"),
                line=dict(
                    # color=smart_format_i["color1"],
                    color=smart_format_i.get("color1", "black"),
                    width=2,
                    )
                )
            )

        return(trace_i)
        #__|

    def __create_layout__(self,
        x_ax_spec="oh",
        title="Scaling Relations",
        showlegend=True,
        ):
        """Create plotly layout dict.

        Args:
            x_ax_spec:
            title:
            showlegend:
        """
        #| - create_layout
        x_ax_title = "G<sub>ads,*OH</sub> (eV)"

        y_ax_title = "G<sub>ads,*OH</sub>, " + \
            "G<sub>ads,*O</sub>, " + \
            "G<sub>ads,*OOH</sub> (eV)"

        tick_lab_size = 12 * (4. / 3.)
        axes_lab_size = 14 * (4. / 3.)
        # legend_size = 18

        #| - Common Axis Dict
        common_axis_dict = {

            # "range": y_axis_range,
            "zeroline": False,
            "showline": True,
            "mirror": 'ticks',
            "linecolor": 'black',
            "showgrid": False,

            "titlefont": dict(size=axes_lab_size),
            "tickfont": dict(
                size=tick_lab_size,
                ),
            "ticks": 'inside',
            "tick0": 0,
            "tickcolor": 'black',
            # "dtick": 0.25,
            "ticklen": 2,
            "tickwidth": 1,
            }
        #__|

        #| - __old__
        # x_range_ooh_vs_oh=[0., 3.5],
        # y_range_ooh_vs_oh=[0., 5.],
        # x_range_o_vs_oh=[0., 3.5],
        # y_range_o_vs_oh=[0., 5.],

        # if y_ax_spec == "ooh":
        #     x_range = self.x_range_ooh_vs_oh
        # elif y_ax_spec == "o":
        #     x_range = self.x_range_o_vs_oh
        # elif y_ax_spec == "oh":
        #     x_range = self.x_range_oh_vs_oh
        # else:
        #     print("Woops - create_layout")
        #
        # if y_ax_spec == "ooh":
        #     y_range = self.y_range_ooh_vs_oh
        # elif y_ax_spec == "o":
        #     y_range = self._range_o_vs_oh
        # elif y_ax_spec == "oh":
        #     y_range = self.y_range_oh_vs_oh
        # else:
        #     print("Woops - create_layout")
        #__|

        x_range = self.x_range
        y_range = self.y_range

        layout_i = {
            "title": title,
            "titlefont": go.layout.Titlefont(size=24),

            "xaxis": dict(
                common_axis_dict,
                **{
                    "title": x_ax_title,
                    "range": x_range,
                    },
                ),

            "yaxis": dict(
                common_axis_dict,
                **{
                    "title": y_ax_title,
                    "range": y_range,
                    },
                ),

            "font": dict(
                family='Arial',
                # size=18,
                color='black',
                ),

            "width": 1.5 * 18.7 * 37.795275591,
            "height": 18.7 * 37.795275591,

            "showlegend": showlegend,

            "legend": dict(
                font=dict(
                    size=10,
                    ),
                ),
            }

        return(layout_i)
        #__|

    def __series_excluded__(self,
        properties_i,
        exclude_dict,
        ):
        """Whether to exclude series_i from fitting.

        Takes an 'exclude_dict' and the series properties_dict and compares
        them key-by-key. If there is a match, then that series is excluded
        (and the function evaluates to True)

        Args:
            properties_i:
            exclude_dict:
        """
        #| - series_excluded
        exclude_dict_keys = list(exclude_dict.keys())
        properties_i_keys = list(properties_i.keys())

        shared_keys = list(
            set(exclude_dict_keys).intersection(set(properties_i_keys)),
            )

        if len(shared_keys) < len(exclude_dict_keys):
            print("series_i doesn't have a specific key!")

        value_match_list = []
        for key_i in shared_keys:
            value_match_i = exclude_dict[key_i] == properties_i[key_i]
            value_match_list.append(value_match_i)


        all_props_match = all(value_match_list)

        # if all_props_match:
        #     print("Ignoring this series for fitting")
        #
        # else:
        #     print("Series not excluded, will include in fitting set")


        return(all_props_match)

        #__|

    def fit_scaling_lines(self,
        dependent_species,
        exclude_dict=None,
        ):
        """Linear fit of either *O or *OOH to *OH

        Args:
            dependent_species:
                y-axis species 'ooh' or 'o'
        """
        #| - fit_scaling_lines

        #| - LOOP
        oh_list = []
        dependent_e_list = []
        for series_i in self.ORR_Free_E_Plot.series_list:

            #| - Excluding series from fitting
            if exclude_dict is not None:
                properties_i = series_i.properties
                exclude_series = self.__series_excluded__(
                    properties_i,
                    exclude_dict,
                    )
                if exclude_series:
                    continue
            #__|

            energy_i = series_i.energy_states_dict[dependent_species]
            dependent_e_list.append(energy_i)
            oh_list.append(series_i.energy_states_dict["oh"])

        #__|

        X = np.array([[i] for i in oh_list])
        y = np.array(dependent_e_list)

        reg = LinearRegression().fit(X, y)

        slope_i = reg.coef_[0]
        intercept_i = reg.intercept_

        out = {"slope": slope_i, "intercept": intercept_i}

        return(out)
        #__|

    def add_ideal_lines(self):
        """Add ideal scaling liknes to plot."""
        #| - add_ideal_lines
        self.add_line({"slope": 1, "intercept": 3.2},
            name="*OOH vs *OH Scaling",
            color="black",
            width=1,
            dash="dash",
            )

        self.add_line({"slope": 2, "intercept": 0.},
            name="*O vs *OH Scaling",
            color="black",
            width=1,
            dash="dash",
            )

        self.add_line({"slope": 1, "intercept": 0.},
            name="*OH vs *OH Scaling",
            color="black",
            width=1,
            dash="dash",
            )
        #__|

    def add_line(self,
        slope_intercept_dict,
        name="add_lines - TEMP",
        color="black",
        width=1,
        dash="dash",
        ):
        """Add line of form y=mx+b to plot.

        Args:
            slope_intercept_dict:
            name:
            color:
            width:
            dash:
        """
        #| - add_line
        slope = slope_intercept_dict["slope"]
        intercept = slope_intercept_dict["intercept"]

        def scaling_meth(E_OH):
            """
            """
            #| - scaling_meth
            out = slope * E_OH + intercept

            return(out)
            #__|

        LH_bound = self.x_range[0]
        RH_bound = self.x_range[1]

        scaling_trace = go.Scatter(
            # x=[self.x_range_ooh_vs_oh[0], self.x_range_ooh_vs_oh[1]],
            x=[LH_bound, RH_bound],
            y=[
                scaling_meth(LH_bound),
                scaling_meth(RH_bound),
                ],
            # name='Fitted scaling',
            name=name,
            mode='lines',
            line=dict(
                dash=dash,
                color=color,
                width=width,
                ),
            )
        # self.data_ooh_oh.append(scaling_trace)
        self.data_lines.append(scaling_trace)
        #__|





    #| - __old__
    # def __ideal_ooh_oh_scaling__(self, E_OH):
    #     """Return the *OOH adsorption energy given DG_*OH by scaling.
    #
    #     Args:
    #         E_OH:DG_*OH energy of adsorption
    #     """
    #     #| - __ideal_ooh_oh_scaling__
    #     return(E_OH + 3.2)
    #     #__|
    #
    # def __ideal_h_oh_scaling__(self, E_OH):
    #     """Return the *OOH adsorption energy given DG_*OH by scaling.
    #
    #     Args:
    #         E_OH: DG_*OH energy of adsorption.
    #     """
    #     #| - __ideal_h_oh_scaling__
    #     return(2 * E_OH)
    #     #__|
    #
    # def __ideal_oh_oh_scaling__(self, E_OH):
    #     """Return the *OH adsorption energy given DG_*OH by scaling.
    #
    #     NOTE: TRIVIAL QUANTITY!!!!!!!!!!!!!!!!!!!
    #
    #     Args:
    #         E_OH: DG_*OH energy of adsorption.
    #     """
    #     #| - __ideal_oh_oh_scaling__
    #     return(E_OH)
    #     #__|
    #
    #__|

    #__| **********************************************************************


class Volcano_Plot():
    """Class to plot OER/ORR volcano plots.

    Development Notes:
        TEMP
    """

    #| - Volcano_Plot *********************************************************

    def __init__(self,
        ORR_Free_E_Plot,
        # x_ax_species="oh",
        x_ax_species="o-oh",  # 'o-oh' or 'oh'
        smart_format_dict=None,

        # mode="all",
        plot_range=None,
        # plot_range={
        #     "y": [1., 5.],
        #     "x": [-2., 4.],
        #     },
        ):
        """
        Input variables to class instance.

        Args:
            ORR_Free_E_Plot:
            mode:
                "all", "ooh_vs_oh", "o_vs_oh"
        """
        #| - __init__
        self.ORR_Free_E_Plot = ORR_Free_E_Plot
        self.x_ax_species = x_ax_species
        self.plot_range = plot_range
        self.smart_format_dict = smart_format_dict

        self.data_points = []
        self.data_lines = []
        #__|

    # NOTE | Rename this create_volcano_plot
    def create_volcano_relations_plot(self,
        smart_format_dict=None,
        ):
        """Create ORR/OER volcano plot.

        Args:
            smart_format_dict:
                Optional dictionary that will format data points
        """
        #| - create_volcano_relations_plot

        #| - Default Smart Format Dict
        smart_format_dict = self.smart_format_dict

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

        #| - Processing Data Points
        x_data_list = []
        y_data_list = []

        for series_i in self.ORR_Free_E_Plot.series_list:

            #| - x-axis energy
            x_spec = self.x_ax_species
            if x_spec == "o-oh":
                e_o = series_i.energy_states_dict["o"]
                e_oh = series_i.energy_states_dict["oh"]
                x_ax_energy = e_o - e_oh
            else:
                x_ax_energy = series_i.energy_states_dict[x_spec]
            #__|

            #| - y-axis limiting potential
            if self.ORR_Free_E_Plot.rxn_type == "ORR":
                lim_pot_i = 1.23 - series_i.overpotential

            elif self.ORR_Free_E_Plot.rxn_type == "OER":
                lim_pot_i = 1.23 + series_i.overpotential_OER
            else:
                print("LSDJFlksdj")
            #__|

            #| - Process series_i
            x_data_list.append(x_ax_energy)
            y_data_list.append(lim_pot_i)

            smart_format_i = self.ORR_Free_E_Plot.__create_smart_format_dict__(
                series_i.properties,
                smart_format_dict,
                )

            # name_i = self.ORR_Free_E_Plot.__create_series_name__(series_i)
            name_i = series_i.series_name
            # print(name_i)
            # print("-------_-_-_-")

            if series_i.color is not None:
                smart_format_i["color2"] = series_i.color

            trace_i = self.__create_trace_i__(
                x_ax_energy,
                lim_pot_i,
                smart_format_i,
                name_i,
                group=series_i.group,
                )

            self.data_points.append(trace_i)
            #__|

        #__|

        #| - Finding plot axis limits
        if self.plot_range is None:
            y_axis_range = [min(y_data_list) - 0.2, max(y_data_list) + 0.2]
            if self.ORR_Free_E_Plot.rxn_type == "OER":
                y_axis_range.reverse()
            else:
                pass

            plot_range = {
                "y": y_axis_range,
                "x": [min(x_data_list) - 0.2, max(x_data_list) + 0.2],
                }

            self.plot_range = plot_range
        #__|

        #__|

    def create_volcano_lines(self,
        gas_molec_dict=None,
        scaling_dict=None,
        plot_all_legs=True,
        plot_min_max_legs=False,
        ):
        """
        """
        #| - create_volcano_lines
        out_data = []

        # x_range = self.plot_range["x"]
        x_range = [-5., 5.]

        #| - Volcano Legs
        volc_legs = [
            'o2_to_ooh',
            'ooh_to_o',
            'o_to_oh',
            'oh_to_h2o',
            ]

        energy_dict = {
            'o2_to_ooh': [],
            'ooh_to_o': [],
            'o_to_oh': [],
            'oh_to_h2o': [],
            }

        x_axis = np.linspace(x_range[0], x_range[1], num=200)
        for leg_i in volc_legs:
            for x_energy_i in x_axis:

                if self.x_ax_species == "oh":
                    g_oh = x_energy_i
                    g_o_minus_g_oh = None

                elif self.x_ax_species == "o-oh":
                    g_oh = None
                    g_o_minus_g_oh = x_energy_i

                energy_dict[leg_i].append(
                    lim_U_i(
                        g_oh=g_oh,
                        g_o_minus_g_oh=g_o_minus_g_oh,
                        # 'o2_to_ooh', 'ooh_to_o', 'o_to_oh', 'oh_to_h2o'
                        mech_step=leg_i,
                        gas_molec_dict=gas_molec_dict,
                        scaling_dict=scaling_dict,
                        rxn_direction="forward",
                        ),
                    )

        if plot_all_legs:
            trace_o2_to_ooh = go.Scatter(
                x=x_axis,
                y=energy_dict["o2_to_ooh"],
                name="O2->*OOH",
                )
            trace_ooh_to_o = go.Scatter(
                x=x_axis,
                y=energy_dict["ooh_to_o"],
                name="*OOH->*O",
                )
            trace_o_to_oh = go.Scatter(
                x=x_axis,
                y=energy_dict["o_to_oh"],
                name="*O->*OH",
                )
            trace_oh_to_h2o = go.Scatter(
                x=x_axis,
                y=energy_dict["oh_to_h2o"],
                name="*OH->H2O",
                )

            out_data.append(trace_o2_to_ooh)
            out_data.append(trace_ooh_to_o)
            out_data.append(trace_o_to_oh)
            out_data.append(trace_oh_to_h2o)
        #__|

        #| - Minimum Energy Legs
        energy_lists = [
            energy_dict["o2_to_ooh"],
            energy_dict["ooh_to_o"],
            energy_dict["o_to_oh"],
            energy_dict["oh_to_h2o"],
            ]

        min_e_list = []
        for leg1, leg2, leg3, leg4 in zip(*energy_lists):
            if self.ORR_Free_E_Plot.rxn_type == "ORR":
                energy_i = min(leg1, leg2, leg3, leg4)

            elif self.ORR_Free_E_Plot.rxn_type == "OER":
                energy_i = max(leg1, leg2, leg3, leg4)

            min_e_list.append(energy_i)

        trace_volcano = go.Scatter(
            x=x_axis,
            y=min_e_list,
            name="activity volcano",
            line=dict(
                color=("black"),
                width=0.5,
                )
            )

        if plot_min_max_legs:
            out_data.append(trace_volcano)

        #__|

        return(out_data)
        #__|

    def __create_trace_i__(self,
        x_energy,
        y_energy,
        smart_format_i,
        name_i,
        # legendgroup=None,
        group=None,
        ):
        """
        """
        #| - __create_trace_i__

        trace_i = go.Scatter(
            x=[x_energy],
            y=[y_energy],
            # mode="markers+text",
            mode="markers",
            name=name_i,
            text=name_i,

            legendgroup=group,

            # textposition='top right',
            textposition='middle left',
            textfont={
                # "family": "Courier New, monospace",
                # "family": font_family,
                "size": 10,
                "color": "black",
                },

            marker=dict(
                size=smart_format_i.get("marker_size", 12),
                color=smart_format_i.get("color2", "red"),
                symbol=smart_format_i.get("symbol", "circle"),
                line=dict(
                    width=smart_format_i.get("marker_border_width", 1.),
                    color=smart_format_i.get("color1", "black"),
                    ),
                ),
            )

        return(trace_i)
        #__|

    def get_plotly_layout(self,
        showlegend=False,
        width=9. * 37.795275591,
        height=9. * 37.795275591,
        ):
        """
        """
        #| - get_plotly_layout

        #| - Properties
        # plot_title="FED"
        plot_title = None
        # plot_title_size = 18
        tick_lab_size = 9 * (4. / 3.)
        axes_lab_size = 10.5 * (4. / 3.)
        legend_size = 18
        # font_family="Computer Modern"  # "Courier New, monospace"
        font_family = "Arial"  # "Courier New, monospace"
        #__|

        # self.x_ax_spec

        if self.x_ax_species == "oh":
            xaxis_title = "dG_*OH (eV)"
        elif self.x_ax_species == "o-oh":
            xaxis_title = "dG_*OH - dG_*O (eV)"

        layout = {
            "title": plot_title,

            "font": {
                "family": font_family,
                "color": "black",
                },

            #| - Axes -----------------------------------------------------

            #| - yaxis
            "yaxis": {
                "title": "Limiting Potential (V)",
                # "title": "$\\Delta G (ev)$",

                "range": self.plot_range["y"],
                "zeroline": False,
                "showline": True,
                "mirror": 'ticks',
                "linecolor": 'black',
                "showgrid": False,

                "titlefont": dict(size=axes_lab_size),

                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                "ticks": 'inside',
                "tick0": 0,
                "tickcolor": 'black',
                "dtick": 0.25,
                "ticklen": 2,
                "tickwidth": 1,
                },
            #__|

            #| - xaxis
            "xaxis": {
                # "title": "$\\Delta G_{OH} (ev)$",
                "title": xaxis_title,
                "range": self.plot_range["x"],
                "zeroline": False,
                "showline": True,
                "mirror": True,
                "linecolor": 'black',
                "showgrid": False,
                "titlefont": dict(size=axes_lab_size),
                "showticklabels": True,
                "ticks": 'inside',
                "tick0": 0,
                "dtick": 0.5,
                "ticklen": 2,
                "tickwidth": 1,
                "tickcolor": 'black',
                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                },
            #__|

            #__|

            # "paper_bgcolor": 'rgba(0,0,0,0)',
            "plot_bgcolor": 'rgba(0,0,0,0)',

            #| - Legend ---------------------------------------------------
            "legend": {
                "traceorder": "normal",
                "font": dict(size=legend_size),
                "x": 0.,
                "y": -0.1,
                # "xanchor": "left",
                "yanchor": "top",
                },

            # "showlegend": False,
            "showlegend": showlegend,

            #__|

            }

        #| - Plot Size Settings
        # bottom_margin_size = 2.5 * 9. * 37.795275591
        plot_size_settings = {
            "width": width,
            "height": height,

            # "width": 9. * 37.795275591,
            # "height": 9 * 37.795275591,

            # "margin": go.layout.Margin({
            #     "l": 50,
            #     "r": 50,
            #     # "b": bottom_margin_size,
            #     # "b": 100,
            #     "b": 1200,
            #     "t": 10,
            #     "pad": 4,
            #     }),
            }

        #__|

        layout = {**layout, **plot_size_settings}

        return(layout)

        #__|


    #__| **********************************************************************


class Free_Energy_Plot():
    """
    Development Notes:
        Take the FED methods out of the ORR_Free_E_Plot class and into this one
    """
    #| - Free_Energy_Plot *****************************************************

    def __init__(self,
        ORR_Free_E_Plot,
        ):
        """
        Input variables to class instance.

        Args:
            ORR_Free_E_Plot:
            mode:
                "all", "ooh_vs_oh", "o_vs_oh"
        """
        #| - __init__
        tmp = 42

        self.ORR_Free_E_Plot = ORR_Free_E_Plot

        # self.x_ax_species = x_ax_species
        # self.smart_format_dict = smart_format_dict
        # self.data_points = []
        # self.data_lines = []
        #__|

    #__| **********************************************************************
