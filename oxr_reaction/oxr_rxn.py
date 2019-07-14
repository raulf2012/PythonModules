#!/usr/bin/env python

"""ORR energetics classes and methods.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression

import plotly.graph_objs as go

pd.options.mode.chained_assignment = None

from oxr_reaction.oxr_series import ORR_Free_E_Series
from oxr_reaction.adsorbate_scaling import lim_U_i
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
        state_title="adsorbate",
        free_e_title="ads_e",
        num_states=5,
        smart_format=None,
        bias=0.,
        color_list=None,
        hover_text_col=None,
        rxn_type="ORR",  # ORR and OER

        # system_properties=None,
        # opt_name=None,
        # properties=None,
        # i_cnt=0,
        # plot_mode="all",
        # smart_format=None,
        # Plotting ************************************************************
        # show_H_e_pairs_annotations=True,
        # show_legend=True,
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

        # self.show_H_e_pairs_annotations = show_H_e_pairs_annotations
        # self.show_legend = show_legend

        # ***********************************
        # COMBAK, moving to FED class, remove later
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

        # self.rxn_x_coord_array = self.create_rxn_coord_array(
        #     self.num_states,
        #     spacing=self.plot_states_sep,
        #     step_size=self.plot_states_width,
        #     )
        #
        # self.mid_state_x_array = self.create_mid_state_x_array()
        # # x_array_data = self.rxn_x_coord_array

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
        format_dict=None,
        overpotential_type="ORR",
        system_properties=None,
        property_key_list=None,
        add_overpot=True,
        color=None,
        overpotential_given=False,
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

        ORR_Series = ORR_Free_E_Series(
            free_energy_df=fe_df,
            properties=system_properties,
            property_key_list=property_key_list,
            state_title=self.state_title,
            free_e_title=self.fe_title,
            group=group,
            bias=self.bias,
            # rxn_x_coord_array=self.rxn_x_coord_array,
            name_i=name_i,
            opt_name=opt_name,  # #######

            color_list=self.color_list,
            color=color,
            hover_text_col=self.hover_text_col,
            plot_mode=plot_mode,  # ##########
            smart_format=smart_format_i,
            format_dict=format_dict,
            add_overpot=add_overpot,
            rxn_type=self.rxn_type,
            overpotential_given=overpotential_given,

            # properties=opt_name,
            # overpotential_type=self.rxn_type,
            # i_cnt=0,  # ##########
            )

        self.series_list.append(ORR_Series)
        #__|

    #__| **********************************************************************






#| - __old__

#
# # ███████ ██████          ██████  ██       ██████  ████████
# # ██      ██   ██         ██   ██ ██      ██    ██    ██
# # ███████ ██████          ██████  ██      ██    ██    ██
# #      ██ ██   ██         ██      ██      ██    ██    ██
# # ███████ ██   ██ ███████ ██      ███████  ██████     ██
#
# class Scaling_Relations_Plot():
#     """Plot scaling relations and some simple fitting schemes.
#
#     Development Notes:
#         IDEA: Add vertical lines to connect *O, *OH, and *OOH data points
#     """
#
#     #| - Scaling_Relations_Plot ***********************************************
#
#     def __init__(self,
#         ORR_Free_E_Plot,
#         mode="all",
#
#         plot_range={
#             "y": [1., 5.],
#             "x": [-2., 4.],
#             },
#
#         x_ax_species="oh",
#
#         marker_color_key="color2",
#         marker_border_color_key="color1",
#         marker_shape_key="symbol",
#         ):
#         """
#         Input variables to class instance.
#
#         Args:
#             ORR_Free_E_Plot:
#             mode:
#                 "all", "ooh_vs_oh", "o_vs_oh"
#         """
#         #| - __init__
#         self.ORR_Free_E_Plot = ORR_Free_E_Plot
#
#         assert (x_ax_species == "oh"), "Only *OH as the x-axis is allowed now"
#         self.x_ax_species = x_ax_species
#         self.marker_color_key = marker_color_key
#         self.marker_border_color_key = marker_border_color_key
#         self.marker_shape_key = marker_shape_key
#
#         # #################################################################
#
#         self.data_points = {
#             "ooh_vs_oh": [],
#             "o_vs_oh": [],
#             "oh_vs_oh": [],
#             }
#         self.data_lines = []
#
#         self.x_range = plot_range["x"]
#         self.y_range = plot_range["y"]
#
#         # self.layout = self.__create_layout__(
#         #     title="Scaling Relations",
#         #     showlegend=True,
#         #     )
#
#         self.scaling_dict = {
#             "ooh": {
#                 "m": None,
#                 "b": None,
#                 },
#
#             "o": {
#                 "m": None,
#                 "b": None,
#                 },
#
#             "oh": {
#                 "m": 1.,
#                 "b": 0.,
#                 },
#
#             }
#
#         self.annotations_list = []
#
#         #__|
#
#     def create_scaling_relations_plot(self,
#         smart_format_dict=None,
#         ):
#         """Return plotly data and layout objects for scaling relations.
#
#         Args:
#             y_ax_spec:
#             x_ax_spec:
#         """
#         #| - create_scaling_relations_plot
#
#         #| - Default Smart Format Dict
#         if smart_format_dict is None:
#             print("No smart format given!")
#             smart_format_dict = [
#                 [{"bulk_system": "IrO3"}, {self.marker_color_key: "green"}],
#                 [{"bulk_system": "IrO2"}, {self.marker_color_key: "yellow"}],
#
#                 [{"coverage_type": "o_covered"}, {self.marker_shape_key: "s"}],
#                 [{"coverage_type": "h_covered"}, {self.marker_shape_key: "^"}],
#
#                 [{"facet": "110"}, {self.marker_border_color_key: "red"}],
#                 [{"facet": "211"}, {self.marker_border_color_key: "green"}],
#                 [{"facet": "100"}, {self.marker_border_color_key: "black"}],
#                 ]
#         #__|
#
#         #| - Processing Data Points
#         for series_i in self.ORR_Free_E_Plot.series_list:
#
#             e_oh = series_i.energy_states_dict["oh"]
#             e_ooh = series_i.energy_states_dict["ooh"]
#             e_o = series_i.energy_states_dict["o"]
#
#
#             # Change self.__create_smart_format_dict__ to,
#             # self.ORR_Free_E_Plot.__create_smart_format_dict__
#             # TODO
#
#             smart_format_i = self.ORR_Free_E_Plot.__create_smart_format_dict__(
#                 series_i.properties,
#                 smart_format_dict,
#                 )
#
#             # This is the old way of getting the name
#             # name_i = self.ORR_Free_E_Plot.__create_series_name__(series_i)
#
#             name_i = series_i.series_name
#
#             if series_i.color is not None:
#                 smart_format_i[self.marker_color_key] = series_i.color
#
#             # NEW, if series has format attached just use that
#             if series_i.format_dict:
#                 smart_format_i = series_i.format_dict
#
#
#             #| - ooh_vs_oh
#             trace_i = self.__create_trace_i__(
#                 e_oh,
#                 e_ooh,
#                 smart_format_i,
#                 name_i,
#                 legendgroup="ooh_vs_oh",
#                 )
#             # self.data_ooh_oh.append(trace_i)
#             self.data_points["ooh_vs_oh"].append(trace_i)
#             #__|
#
#             #| - o_vs_oh
#             trace_i = self.__create_trace_i__(
#                 e_oh,
#                 e_o,
#                 smart_format_i,
#                 name_i,
#                 legendgroup="o_vs_oh",
#                 )
#             # self.data_o_oh.append(trace_i)
#             self.data_points["o_vs_oh"].append(trace_i)
#             #__|
#
#             #| - oh_vs_oh
#             trace_i = self.__create_trace_i__(
#                 e_oh,
#                 e_oh,
#                 smart_format_i,
#                 name_i,
#                 legendgroup="oh_vs_oh",
#                 )
#             # self.data_oh_oh.append(trace_i)
#             self.data_points["oh_vs_oh"].append(trace_i)
#             #__|
#
#         #__|
#
#         #__|
#
#     # Deprecated, delete this later
#     def __create_smart_format_dict__(self, property_dict, smart_format_dict):
#         """Create smart format dictionary.
#
#         Args:
#             property_dict:
#             smart_format_dict:
#         """
#         #| - __create_smart_format_dict__
#         format_dict = {}
#         for key_i, value_i in property_dict.items():
#             for format_i in smart_format_dict:
#                 if list(format_i[0])[0] == key_i:
#                     if list(format_i[0].values())[0] == value_i:
#                         format_dict.update(format_i[1])
#
#         return(format_dict)
#         #__|
#
#     def __create_series_name__(self, series_i):
#         """
#         """
#         #| - create_series_name
#         name_i = ""
#         for key, value in series_i.properties.items():
#             if key == "coverage":
#                 continue
#
#             name_i += str(key) + ": " + str(value) + " | "
#
#         return(name_i)
#         #__|
#
#     def __create_trace_i__(self,
#         x_energy,
#         y_energy,
#         smart_format_i,
#         name_i,
#         legendgroup=None,
#         ):
#         """
#         """
#         #| - create_trace_i
#         # NOTE Looks like I need to put these in a list here
#         x_energy = [x_energy]
#         y_energy = [y_energy]
#
#         trace_i = go.Scatter(
#             x=x_energy,
#             y=y_energy,
#             text=name_i,
#             name=name_i,
#             mode="markers",
#             legendgroup=legendgroup,
#             marker=dict(
#                 size=smart_format_i.get("marker_size", 9),
#                 symbol=smart_format_i.get(
#                     self.marker_shape_key, "circle"),
#                 color=smart_format_i.get(
#                     self.marker_color_key, "pink"),
#                 line=dict(
#                     # color=smart_format_i[marker_border_color_key],
#                     color=smart_format_i.get(
#                         self.marker_border_color_key, "black"),
#                     width=1.,
#                     )
#                 )
#             )
#
#         return(trace_i)
#         #__|
#
#     # NOTE | This shouldn't be an internal method
#     def __create_layout__(self,
#         # x_ax_spec="oh",
#         title="Scaling Relations",
#         showlegend=True,
#         layout_dict=None,
#         ):
#         """Create plotly layout dict.
#
#         Args:
#             x_ax_spec:
#             title:
#             showlegend:
#         """
#         #| - create_layout
#
#         # if x_ax_spec == ""
#         if self.x_ax_species == "oh":
#             x_ax_title = "G<sub>ads,*OH</sub> (eV)"
#         else:
#             x_ax_title = "TEMP"
#
#         y_ax_title = "G<sub>ads,*OH</sub>, " + \
#             "G<sub>ads,*O</sub>, " + \
#             "G<sub>ads,*OOH</sub> (eV)"
#
#         tick_lab_size = 12 * (4. / 3.)
#         axes_lab_size = 14 * (4. / 3.)
#
#         # legend_size = 18
#
#         #| - Common Axis Dict
#         common_axis_dict = {
#
#             # "range": y_axis_range,
#             "zeroline": False,
#             "showline": True,
#             "mirror": 'ticks',
#             "linecolor": 'black',
#             "showgrid": False,
#
#             "titlefont": dict(size=axes_lab_size),
#             "tickfont": dict(
#                 size=tick_lab_size,
#                 ),
#             "ticks": 'inside',
#             "tick0": 0,
#             "tickcolor": 'black',
#             # "dtick": 0.25,
#             "ticklen": 2,
#             "tickwidth": 1,
#             }
#         #__|
#
#         #| - __old__
#         # x_range_ooh_vs_oh=[0., 3.5],
#         # y_range_ooh_vs_oh=[0., 5.],
#         # x_range_o_vs_oh=[0., 3.5],
#         # y_range_o_vs_oh=[0., 5.],
#
#         # if y_ax_spec == "ooh":
#         #     x_range = self.x_range_ooh_vs_oh
#         # elif y_ax_spec == "o":
#         #     x_range = self.x_range_o_vs_oh
#         # elif y_ax_spec == "oh":
#         #     x_range = self.x_range_oh_vs_oh
#         # else:
#         #     print("Woops - create_layout")
#         #
#         # if y_ax_spec == "ooh":
#         #     y_range = self.y_range_ooh_vs_oh
#         # elif y_ax_spec == "o":
#         #     y_range = self._range_o_vs_oh
#         # elif y_ax_spec == "oh":
#         #     y_range = self.y_range_oh_vs_oh
#         # else:
#         #     print("Woops - create_layout")
#         #__|
#
#         x_range = self.x_range
#         y_range = self.y_range
#
#         layout_i = {
#             "title": title,
#             "titlefont": go.layout.title.Font(size=24),
#             # "titlefont": go.layout.Titlefont(size=24),
#
#             "xaxis": dict(
#                 common_axis_dict,
#                 **{
#                     "title": x_ax_title,
#                     "range": x_range,
#                     },
#                 ),
#
#             "yaxis": dict(
#                 common_axis_dict,
#                 **{
#                     "title": y_ax_title,
#                     "range": y_range,
#                     },
#                 ),
#
#             # Margin
#             "margin": go.layout.Margin(
#                 b=50.,
#                 l=50.,
#                 r=50.,
#                 t=50.,
#                 ),
#
#             "font": dict(
#                 family='Arial',
#                 # size=18,
#                 color='black',
#                 ),
#
#             "width": 1.5 * 18.7 * 37.795275591,
#             "height": 18.7 * 37.795275591,
#
#             "showlegend": showlegend,
#
#             "legend": dict(
#                 font=dict(
#                     size=10,
#                     ),
#                 ),
#             }
#
#         if layout_dict is not None:
#             from misc_modules.misc_methods import dict_merge
#             dict_merge(layout_i, layout_dict)
#             # layout_i = {**layout_i, **layout_dict}
#
#         return(layout_i)
#         #__|
#
#     def __series_excluded__(self,
#         properties_i,
#         exclude_dict,
#         ):
#         """Whether to exclude series_i from fitting.
#
#         Takes an 'exclude_dict' and the series properties_dict and compares
#         them key-by-key. If there is a match, then that series is excluded
#         (and the function evaluates to True)
#
#         Args:
#             properties_i:
#             exclude_dict:
#         """
#         #| - series_excluded
#         exclude_dict_keys = list(exclude_dict.keys())
#         properties_i_keys = list(properties_i.keys())
#
#         shared_keys = list(
#             set(exclude_dict_keys).intersection(set(properties_i_keys)),
#             )
#
#         if len(shared_keys) < len(exclude_dict_keys):
#             print("series_i doesn't have a specific key!")
#
#         value_match_list = []
#         for key_i in shared_keys:
#             value_match_i = exclude_dict[key_i] == properties_i[key_i]
#             value_match_list.append(value_match_i)
#
#
#         all_props_match = all(value_match_list)
#
#         # if all_props_match:
#         #     print("Ignoring this series for fitting")
#         #
#         # else:
#         #     print("Series not excluded, will include in fitting set")
#
#
#         return(all_props_match)
#
#         #__|
#
#     def fit_scaling_lines(self,
#         dependent_species,  # 'ooh', 'o', 'oh'
#         exclude_dict=None,
#         ):
#         """Linear fit of either *O or *OOH to *OH
#
#         Args:
#             dependent_species:
#                 y-axis species 'ooh' or 'o'
#         """
#         #| - fit_scaling_lines
#
#         #| - LOOP
#         oh_list = []
#         dependent_e_list = []
#         for series_i in self.ORR_Free_E_Plot.series_list:
#
#             #| - Excluding series from fitting
#             if exclude_dict is not None:
#                 properties_i = series_i.properties
#                 exclude_series = self.__series_excluded__(
#                     properties_i,
#                     exclude_dict,
#                     )
#                 if exclude_series:
#                     continue
#             #__|
#
#             energy_i = series_i.energy_states_dict[dependent_species]
#             dependent_e_list.append(energy_i)
#             oh_list.append(series_i.energy_states_dict["oh"])
#
#         #__|
#
#         X = np.array([[i] for i in oh_list])
#         y = np.array(dependent_e_list)
#
#         reg = LinearRegression().fit(X, y)
#
#         slope_i = reg.coef_[0]
#         intercept_i = reg.intercept_
#
#         print("Scaling fit for ", dependent_species)
#         print("intercept_i: ", str(intercept_i))
#         print("slope_i: ", str(slope_i))
#         print("")
#
#         out = {"slope": slope_i, "intercept": intercept_i}
#
#         self.scaling_dict[dependent_species] = {
#             "m": slope_i,
#             "b": intercept_i,
#             }
#         # print("_------__)_Z(*XF(8))")
#
#         #| - Equation Annotations
#         if dependent_species == "ooh":
#             eqn_str_i = ("" +
#                 "G<sub>OOH</sub>=" +
#                 str(round(slope_i, 4)) +
#                 " G<sub>OH</sub>+" +
#                 str(round(intercept_i, 4)) +
#                 ""
#                 )
#
#         elif dependent_species == "o":
#             eqn_str_i = ("" +
#                 "G<sub>O</sub> = " +
#                 str(round(slope_i, 4)) +
#                 " G<sub>OH</sub>+" +
#                 str(round(intercept_i, 4)) +
#                 ""
#                 )
#
#         elif dependent_species == "oh":
#             eqn_str_i = ("" +
#                 "G<sub>OH</sub> = " +
#                 str(round(slope_i, 4)) +
#                 " G<sub>OH</sub>+" +
#                 str(round(intercept_i, 4)) +
#                 ""
#                 )
#
#         else:
#             eqn_str_i = "TEMP TEMP TEMP TEMP | 190213 | RF"
#             raise ValueError('A very specific bad thing happened.')
#
#         annotation_i = dict(
#             x=0.,
#             y=1.,
#             xref="paper",
#             yref="paper",
#             text=eqn_str_i,
#             font=dict(
#                 color="red",
#                 family="Droid Sans Mono,Overpass",
#                 size=9. * (4. / 3.),
#                 ),
#             showarrow=False,
#             xanchor="left",
#             yshift=-11. * (4. / 3.) * len(self.annotations_list),
#             xshift=5.,
#             )
#
#         self.annotations_list.append(annotation_i)
#         #__|
#
#
#         return(out)
#         #__|
#
#     def add_ideal_lines(self):
#         """Add ideal scaling liknes to plot."""
#         #| - add_ideal_lines
#         self.add_line({"slope": 1, "intercept": 3.2},
#             name="*OOH vs *OH Scaling",
#             color="black",
#             width=1,
#             dash="dash",
#             )
#
#         self.add_line({"slope": 2, "intercept": 0.},
#             name="*O vs *OH Scaling",
#             color="black",
#             width=1,
#             dash="dash",
#             )
#
#         self.add_line({"slope": 1, "intercept": 0.},
#             name="*OH vs *OH Scaling",
#             color="black",
#             width=1,
#             dash="dash",
#             )
#         #__|
#
#     def add_line(self,
#         slope_intercept_dict,
#         name="add_lines - TEMP",
#         color="black",
#         width=1,
#         dash="dash",
#         ):
#         """Add line of form y=mx+b to plot.
#
#         Args:
#             slope_intercept_dict:
#             name:
#             color:
#             width:
#             dash:
#         """
#         #| - add_line
#
#         # print(slope_intercept_dict)
#
#         slope = slope_intercept_dict["slope"]
#         intercept = slope_intercept_dict["intercept"]
#
#         def scaling_meth(E_OH):
#             """
#             """
#             #| - scaling_meth
#             out = slope * E_OH + intercept
#
#             return(out)
#             #__|
#
#         LH_bound = self.x_range[0]
#         RH_bound = self.x_range[1]
#
#         scaling_trace = go.Scatter(
#             # x=[self.x_range_ooh_vs_oh[0], self.x_range_ooh_vs_oh[1]],
#             x=[LH_bound, RH_bound],
#             y=[
#                 scaling_meth(LH_bound),
#                 scaling_meth(RH_bound),
#                 ],
#             # name='Fitted scaling',
#             name=name,
#             mode='lines',
#             line=dict(
#                 dash=dash,
#                 color=color,
#                 width=width,
#                 ),
#             )
#         # self.data_ooh_oh.append(scaling_trace)
#         self.data_lines.append(scaling_trace)
#
#         #
#         # # Annotation
#         # ooh_vs_oh_eqn = ("" +
#         #     "G_*OOH = " +
#         #     str(round(SC_PLT.scaling_dict["ooh"]["m"], 5)) +
#         #     " G_*OH + " +
#         #     str(round(SC_PLT.scaling_dict["ooh"]["b"], 5)) +
#         #     ""
#         #     )
#         #
#         # o_vs_oh_eqn = ("" +
#         #     "G_*O  =  " +
#         #     str(round(SC_PLT.scaling_dict["o"]["m"], 5)) +
#         #     " G_*OH + " +
#         #     str(round(SC_PLT.scaling_dict["o"]["b"], 5)) +
#         #     ""
#         #     )
#
#         #__|
#
#
#
#
#
#     #| - __old__
#     # def __ideal_ooh_oh_scaling__(self, E_OH):
#     #     """Return the *OOH adsorption energy given DG_*OH by scaling.
#     #
#     #     Args:
#     #         E_OH:DG_*OH energy of adsorption
#     #     """
#     #     #| - __ideal_ooh_oh_scaling__
#     #     return(E_OH + 3.2)
#     #     #__|
#     #
#     # def __ideal_h_oh_scaling__(self, E_OH):
#     #     """Return the *OOH adsorption energy given DG_*OH by scaling.
#     #
#     #     Args:
#     #         E_OH: DG_*OH energy of adsorption.
#     #     """
#     #     #| - __ideal_h_oh_scaling__
#     #     return(2 * E_OH)
#     #     #__|
#     #
#     # def __ideal_oh_oh_scaling__(self, E_OH):
#     #     """Return the *OH adsorption energy given DG_*OH by scaling.
#     #
#     #     NOTE: TRIVIAL QUANTITY!!!!!!!!!!!!!!!!!!!
#     #
#     #     Args:
#     #         E_OH: DG_*OH energy of adsorption.
#     #     """
#     #     #| - __ideal_oh_oh_scaling__
#     #     return(E_OH)
#     #     #__|
#     #
#     #__|
#
# #__| **********************************************************************




# # ██    ██  ██████  ██       ██████         ██████  ██       ██████  ████████
# # ██    ██ ██    ██ ██      ██              ██   ██ ██      ██    ██    ██
# # ██    ██ ██    ██ ██      ██              ██████  ██      ██    ██    ██
# #  ██  ██  ██    ██ ██      ██              ██      ██      ██    ██    ██
# #   ████    ██████  ███████  ██████ ███████ ██      ███████  ██████     ██
#
# class Volcano_Plot():
#     """Class to plot OER/ORR volcano plots.
#
#     Development Notes:
#         TEMP
#     """
#
#     #| - Volcano_Plot *********************************************************
#
#     def __init__(self,
#         ORR_Free_E_Plot,
#         x_ax_species="o-oh",  # 'o-oh' or 'oh'
#         smart_format_dict=None,
#         plot_range=None,
#         marker_color_key="color2",
#         marker_border_color_key="color1",
#         marker_shape_key="symbol",
#         ):
#         """
#         Input variables to class instance.
#
#         Args:
#             ORR_Free_E_Plot:
#             mode:
#                 "all", "ooh_vs_oh", "o_vs_oh"
#             plot_range:
#                 Ex.)
#                 plot_range={
#                     "y": [1., 5.],
#                     "x": [-2., 4.],
#                     }
#
#         """
#         #| - __init__
#         self.ORR_Free_E_Plot = ORR_Free_E_Plot
#         self.x_ax_species = x_ax_species
#         self.plot_range = plot_range
#         self.smart_format_dict = smart_format_dict
#
#         self.data_points = []
#         self.data_lines = []
#
#         self.marker_color_key = marker_color_key
#         self.marker_border_color_key = marker_border_color_key
#         self.marker_shape_key = marker_shape_key
#
#         #__|
#
#     # NOTE | Rename this create_volcano_plot
#     def create_volcano_relations_plot(self,
#         show_data_labels=False,
#         # smart_format_dict=None,
#         ):
#         """Create ORR/OER volcano plot.
#
#         Args:
#             smart_format_dict:
#                 Optional dictionary that will format data points
#         """
#         #| - create_volcano_relations_plot
#
#         #| - Default Smart Format Dict
#         smart_format_dict = self.smart_format_dict
#
#         if smart_format_dict is None:
#             print("No smart format given!")
#             smart_format_dict = [
#                 [{"bulk_system": "IrO3"}, {self.marker_color_key: "green"}],
#                 [{"bulk_system": "IrO2"}, {self.marker_color_key: "yellow"}],
#
#                 [{"coverage_type": "o_covered"}, {"symbol": "s"}],
#                 [{"coverage_type": "h_covered"}, {"symbol": "^"}],
#
#                 [{"facet": "110"}, {"color1": "red"}],
#                 [{"facet": "211"}, {"color1": "green"}],
#                 [{"facet": "100"}, {"color1": "black"}],
#                 ]
#         #__|
#
#         #| - Processing Data Points
#         x_data_list = []
#         y_data_list = []
#
#         for series_i in self.ORR_Free_E_Plot.series_list:
#
#             #| - x-axis energy
#             x_spec = self.x_ax_species
#             if x_spec == "o-oh":
#                 e_o = series_i.energy_states_dict["o"]
#                 e_oh = series_i.energy_states_dict["oh"]
#                 x_ax_energy = e_o - e_oh
#             else:
#                 x_ax_energy = series_i.energy_states_dict[x_spec]
#             #__|
#
#             #| - y-axis limiting potential
#             if self.ORR_Free_E_Plot.rxn_type == "ORR":
#                 lim_pot_i = 1.23 - series_i.overpotential
#
#             elif self.ORR_Free_E_Plot.rxn_type == "OER":
#                 lim_pot_i = 1.23 + series_i.overpotential_OER
#             else:
#                 print("LSDJFlksdj")
#             #__|
#
#             #| - Process series_i
#             x_data_list.append(x_ax_energy)
#             y_data_list.append(lim_pot_i)
#
#             smart_format_i = self.ORR_Free_E_Plot.__create_smart_format_dict__(
#                 series_i.properties,
#                 smart_format_dict,
#                 )
#
#             name_i = series_i.series_name
#
#             if series_i.color is not None:
#                 smart_format_i[self.marker_color_key] = series_i.color
#
#
#             format_i = smart_format_i
#
#             if series_i.format_dict:
#                 format_i = series_i.format_dict
#
#             trace_i = self.__create_trace_i__(
#                 x_ax_energy,
#                 lim_pot_i,
#                 # smart_format_i,
#                 format_i,
#                 name_i,
#                 group=series_i.group,
#                 show_data_labels=show_data_labels,
#                 )
#
#             self.data_points.append(trace_i)
#             #__|
#
#         #__|
#
#         #| - Finding plot axis limits
#         if self.plot_range is None:
#             y_axis_range = [min(y_data_list) - 0.2, max(y_data_list) + 0.2]
#             if self.ORR_Free_E_Plot.rxn_type == "OER":
#                 y_axis_range.reverse()
#             else:
#                 pass
#
#             plot_range = {
#                 "y": y_axis_range,
#                 "x": [min(x_data_list) - 0.2, max(x_data_list) + 0.2],
#                 }
#
#             self.plot_range = plot_range
#         #__|
#
#         #__|
#
#     def create_volcano_lines(self,
#         gas_molec_dict=None,
#         scaling_dict=None,
#         plot_all_legs=True,
#         plot_min_max_legs=False,
#         trace_priority="top",  # 'top' or 'bottom'
#         legs_to_plot=[
#             "o2_to_ooh",
#             "ooh_to_o",
#             "o_to_oh",
#             "oh_to_h2o",
#             ],
#         line_color="black",
#         ):
#         """Create volcano data traces.
#
#         Args:
#             gas_molec_dict:
#             scaling_dict:
#             plot_all_legs:
#             plot_min_max_legs:
#             trace_priority:
#                 if 'top', the volcano lines will be placed on the top of the
#                 plot, if 'bottom' then the data points will by on top of the
#                 volcano
#         """
#         #| - create_volcano_lines
#         out_data = []
#         x_range = self.plot_range["x"]
#
#         #| - Volcano Legs
#         volc_legs = [
#             'o2_to_ooh',
#             'ooh_to_o',
#             'o_to_oh',
#             'oh_to_h2o',
#             ]
#
#         energy_dict = {
#             'o2_to_ooh': [],
#             'ooh_to_o': [],
#             'o_to_oh': [],
#             'oh_to_h2o': [],
#             }
#
#         #| - Create Volcano Legs (LOOP)
#         x_axis = np.linspace(x_range[0], x_range[1], num=500)
#         for leg_i in volc_legs:
#             for x_energy_i in x_axis:
#
#                 if self.x_ax_species == "oh":
#                     g_oh = x_energy_i
#                     g_o_minus_g_oh = None
#
#                 elif self.x_ax_species == "o-oh":
#                     g_oh = None
#                     g_o_minus_g_oh = x_energy_i
#
#                 energy_dict[leg_i].append(
#                     lim_U_i(
#                         g_oh=g_oh,
#                         g_o_minus_g_oh=g_o_minus_g_oh,
#                         # 'o2_to_ooh', 'ooh_to_o', 'o_to_oh', 'oh_to_h2o'
#                         mech_step=leg_i,
#                         gas_molec_dict=gas_molec_dict,
#                         scaling_dict=scaling_dict,
#                         rxn_direction="forward",
#                         ),
#                     )
#         #__|
#
#         if plot_all_legs:
#
#             #| - plot_all_legs
#             # hoverinfo_type = "none"
#             hoverinfo_type = "name"
#
#             trace_o2_to_ooh = go.Scatter(
#                 x=x_axis,
#                 y=energy_dict["o2_to_ooh"],
#                 name="O2->*OOH",
#                 hoverinfo=hoverinfo_type,
#                 line=dict(
#                     color="#e7b8bc",
#                     width=2,
#                     dash="solid",
#                     )
#                 )
#
#             trace_ooh_to_o = go.Scatter(
#                 x=x_axis,
#                 y=energy_dict["ooh_to_o"],
#                 name="*OOH->*O",
#                 hoverinfo=hoverinfo_type,
#                 line=dict(
#                     color="#afd7c3",
#                     width=2,
#                     dash="solid",
#                     )
#                 )
#
#             trace_o_to_oh = go.Scatter(
#                 x=x_axis,
#                 y=energy_dict["o_to_oh"],
#                 name="*O->*OH",
#                 hoverinfo=hoverinfo_type,
#                 line=dict(
#                     color="#b5c4e2",
#                     width=2,
#                     dash="solid",
#                     )
#                 )
#
#             trace_oh_to_h2o = go.Scatter(
#                 x=x_axis,
#                 y=energy_dict["oh_to_h2o"],
#                 name="*OH->H2O",
#                 hoverinfo=hoverinfo_type,
#                 line=dict(
#                     color="#dbcdab",
#                     width=2,
#                     dash="solid",
#                     )
#                 )
#
#             if trace_priority == "top":
#                 out_data.append(trace_o2_to_ooh)
#                 out_data.append(trace_ooh_to_o)
#                 out_data.append(trace_o_to_oh)
#                 out_data.append(trace_oh_to_h2o)
#
#             elif trace_priority == "bottom":
#                 out_data.insert(0, trace_o2_to_ooh)
#                 out_data.insert(0, trace_ooh_to_o)
#                 out_data.insert(0, trace_o_to_oh)
#                 out_data.insert(0, trace_oh_to_h2o)
#             #__|
#
#         #__|
#
#         #| - Minimum Energy Legs
#         energy_lists= []
#         for leg_i in legs_to_plot:
#             energy_lists.append(energy_dict[leg_i])
#
#         min_max_e_list = []
#         for legs in zip(*energy_lists):
#             if self.ORR_Free_E_Plot.rxn_type == "ORR":
#                 energy_i = min(*legs)
#
#             elif self.ORR_Free_E_Plot.rxn_type == "OER":
#                 energy_i = max(*legs)
#
#             min_max_e_list.append(energy_i)
#
#         trace_volcano = go.Scatter(
#             x=x_axis,
#             y=min_max_e_list,
#             name="activity volcano",
#             hoverinfo="skip",
#             line=dict(
#                 color=line_color,
#                 width=2,
#                 # dash="dash",
#                 dash="5px,2px,5px,2px",
#                 )
#             )
#
#         if plot_min_max_legs:
#             if trace_priority == "top":
#                 out_data.append(trace_volcano)
#
#             elif trace_priority == "bottom":
#                 out_data.insert(0, trace_volcano)
#         #__|
#
#         return(out_data)
#         #__|
#
#     def __create_trace_i__(self,
#         x_energy,
#         y_energy,
#         smart_format_i,
#         name_i,
#         # legendgroup=None,
#         group=None,
#         show_data_labels=False,
#         ):
#         """
#         """
#         #| - __create_trace_i__
#
#         if show_data_labels is True:
#             mode_i = "markers+text"
#         elif show_data_labels is False:
#             mode_i = "markers"
#         else:
#             print("TEMP TEMP TEMP | __create_trace_i__")
#
#         # print(mode_i)
#
#         trace_i = go.Scatter(
#             x=[x_energy],
#             y=[y_energy],
#
#             mode=mode_i,
#             # mode="markers+text",
#             # mode="markers",
#
#             name=name_i,
#             text=[name_i],
#             # text=["TEMP"],
#
#             legendgroup=group,
#
#             hoverinfo="text",
#
#             # hoverinfo=None,
#             # hoverinfosrc=None,
#             # hoverlabel=None,
#             # hoveron=None,
#             # hovertext=None,
#             # hovertextsrc=None,
#
#             # textposition='top right',
#             textposition='middle left',
#             textfont={
#                 # "family": "Courier New, monospace",
#                 # "family": font_family,
#                 "size": 10,
#                 "color": "black",
#                 },
#
#             marker=dict(
#                 size=smart_format_i.get("marker_size", 9),
#                 color=smart_format_i.get(self.marker_color_key, "red"),
#                 symbol=smart_format_i.get(
#                     self.marker_shape_key, "circle"),
#                 line=dict(
#                     width=smart_format_i.get("marker_border_width", 1.),
#                     color=smart_format_i.get(
#                         self.marker_border_color_key, "black"),
#                     ),
#                 ),
#             )
#
#         return(trace_i)
#         #__|
#
#     def get_plotly_layout(self,
#         showlegend=False,
#         width=9. * 37.795275591,
#         height=9. * 37.795275591,
#         layout_dict=None,
#         ):
#         """
#         """
#         #| - get_plotly_layout
#
#         #| - Properties
#         # plot_title="FED"
#         plot_title = None
#         # plot_title_size = 18
#         # tick_lab_size = 9 * (4. / 3.)
#         tick_lab_size = 8 * (4. / 3.)
#         axes_lab_size = 9 * (4. / 3.)
#         legend_size = 18
#         # font_family="Computer Modern"  # "Courier New, monospace"
#         font_family = "Arial"  # "Courier New, monospace"
#         #__|
#
#         # self.x_ax_spec
#
#         if self.x_ax_species == "oh":
#             # xaxis_title = "dG_*OH (eV)"
#             xaxis_title = "dG<sub>OH</sub> (eV)"
#         elif self.x_ax_species == "o-oh":
#             # xaxis_title = "dG_*OH - dG_*O (eV)"
#             xaxis_title = "dG<sub>O</sub> - dG<sub>OH</sub> (eV)"
#
#         # layout["margin"] = go.layout.Margin(
#         #     b=50.,
#         #     l=50.,
#         #     r=50.,
#         #     t=50.,
#         # #     pad=20.,
#         #     )
#
#         layout = {
#             "title": plot_title,
#
#             "font": {
#                 "family": font_family,
#                 "color": "black",
#                 },
#
#             #| - Axes -----------------------------------------------------
#
#             #| - yaxis
#             "yaxis": {
#                 "title": "Limiting Potential (V)",
#                 # "title": "$\\Delta G (ev)$",
#
#                 "range": self.plot_range["y"],
#                 "zeroline": False,
#                 "showline": True,
#                 "mirror": 'ticks',
#                 # "linecolor": 'red',
#                 "linecolor": 'black',
#                 "showgrid": False,
#
#                 "titlefont": dict(size=axes_lab_size),
#
#                 "tickfont": dict(
#                     size=tick_lab_size,
#                     ),
#                 "ticks": 'inside',
#                 "tick0": 0,
#                 "tickcolor": 'black',
#                 "dtick": 0.1,
#                 "ticklen": 2,
#                 "tickwidth": 1,
#                 },
#             #__|
#
#             #| - xaxis
#             "xaxis": {
#                 # "title": "$\\Delta G_{OH} (ev)$",
#                 "title": xaxis_title,
#                 "range": self.plot_range["x"],
#                 "zeroline": False,
#                 "showline": True,
#                 "mirror": True,
#                 # "linecolor": 'red',
#                 "linecolor": 'black',
#                 "showgrid": False,
#                 "titlefont": dict(size=axes_lab_size),
#                 "showticklabels": True,
#                 "ticks": 'inside',
#                 "tick0": 0,
#                 "dtick": 0.2,
#                 "ticklen": 2,
#                 "tickwidth": 1,
#                 "tickcolor": 'black',
#                 "tickfont": dict(
#                     size=tick_lab_size,
#                     ),
#                 },
#             #__|
#
#             #__|
#
#             "margin": go.layout.Margin(
#                 b=50.,
#                 l=50.,
#                 r=50.,
#                 t=50.,
#                 ),
#
#             # "paper_bgcolor": 'rgba(0,0,0,0)',
#             "plot_bgcolor": 'rgba(0,0,0,0)',
#
#             #| - Legend ---------------------------------------------------
#             "legend": {
#                 "traceorder": "normal",
#                 "font": dict(size=legend_size),
#                 "x": 0.,
#                 "y": -0.1,
#                 # "xanchor": "left",
#                 "yanchor": "top",
#                 },
#
#             # "showlegend": False,
#             "showlegend": showlegend,
#
#             #__|
#
#             }
#
#         #| - Plot Size Settings
#         # bottom_margin_size = 2.5 * 9. * 37.795275591
#         plot_size_settings = {
#             "width": width,
#             "height": height,
#
#             # "width": 9. * 37.795275591,
#             # "height": 9 * 37.795275591,
#
#             # "margin": go.layout.Margin({
#             #     "l": 50,
#             #     "r": 50,
#             #     # "b": bottom_margin_size,
#             #     # "b": 100,
#             #     "b": 1200,
#             #     "t": 10,
#             #     "pad": 4,
#             #     }),
#             }
#
#         #__|
#
#         layout = {**layout, **plot_size_settings}
#
#         #| - Applying Layout override dict
#         if layout_dict is not None:
#             from misc_modules.misc_methods import dict_merge
#             dict_merge(layout, layout_dict)
#
#             # layout_i = {**layout_i, **layout_dict}
#
#         #__|
#
#         return(layout)
#
#         #__|
#
#
#     #__| **********************************************************************
#__|
