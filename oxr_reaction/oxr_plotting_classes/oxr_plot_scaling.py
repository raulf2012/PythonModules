#!/usr/bin/env python

"""OXR class for plotting OXR energetics (scaling relations).

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression

import plotly.graph_objs as go

# pd.options.mode.chained_assignment = None

from oxr_reaction.oxr_series import ORR_Free_E_Series
from oxr_reaction.adsorbate_scaling import lim_U_i
#__|


# ███████ ██████          ██████  ██       ██████  ████████
# ██      ██   ██         ██   ██ ██      ██    ██    ██
# ███████ ██████          ██████  ██      ██    ██    ██
#      ██ ██   ██         ██      ██      ██    ██    ██
# ███████ ██   ██ ███████ ██      ███████  ██████     ██

class Scaling_Relations_Plot():
    """Plot scaling relations and some simple fitting schemes.

    Development Notes:
        IDEA: Add vertical lines to connect *O, *OH, and *OOH data points
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

        marker_color_key="color2",
        marker_border_color_key="color1",
        marker_shape_key="symbol",
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
        self.x_ax_species = x_ax_species
        self.marker_color_key = marker_color_key
        self.marker_border_color_key = marker_border_color_key
        self.marker_shape_key = marker_shape_key

        # #################################################################

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

        self.scaling_dict = {
            "ooh": {
                "m": None,
                "b": None,
                },

            "o": {
                "m": None,
                "b": None,
                },

            "oh": {
                "m": 1.,
                "b": 0.,
                },

            }

        self.annotations_list = []

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
                [{"bulk_system": "IrO3"}, {self.marker_color_key: "green"}],
                [{"bulk_system": "IrO2"}, {self.marker_color_key: "yellow"}],

                [{"coverage_type": "o_covered"}, {self.marker_shape_key: "s"}],
                [{"coverage_type": "h_covered"}, {self.marker_shape_key: "^"}],

                [{"facet": "110"}, {self.marker_border_color_key: "red"}],
                [{"facet": "211"}, {self.marker_border_color_key: "green"}],
                [{"facet": "100"}, {self.marker_border_color_key: "black"}],
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
                smart_format_i[self.marker_color_key] = series_i.color

            # NEW, if series has format attached just use that
            if series_i.format_dict:
                smart_format_i = series_i.format_dict


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
            # legendgroup=legendgroup,
            legendgroup=name_i,
            marker=dict(
                size=smart_format_i.get("marker_size", 9),
                symbol=smart_format_i.get(
                    self.marker_shape_key, "circle"),
                color=smart_format_i.get(
                    self.marker_color_key, "pink"),
                line=dict(
                    # color=smart_format_i[marker_border_color_key],
                    color=smart_format_i.get(
                        self.marker_border_color_key, "black"),
                    width=1.,
                    )
                )
            )

        return(trace_i)
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
        dependent_species,  # 'ooh', 'o', 'oh'
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

        print("Scaling fit for ", dependent_species)
        print("intercept_i: ", str(intercept_i))
        print("slope_i: ", str(slope_i))
        print("")

        out = {"slope": slope_i, "intercept": intercept_i}

        self.scaling_dict[dependent_species] = {
            "m": slope_i,
            "b": intercept_i,
            }
        # print("_------__)_Z(*XF(8))")

        #| - Equation Annotations
        if dependent_species == "ooh":
            eqn_str_i = ("" +
                "G<sub>OOH</sub>=" +
                str(round(slope_i, 4)) +
                " G<sub>OH</sub>+" +
                str(round(intercept_i, 4)) +
                ""
                )

        elif dependent_species == "o":
            eqn_str_i = ("" +
                "G<sub>O</sub> = " +
                str(round(slope_i, 4)) +
                " G<sub>OH</sub>+" +
                str(round(intercept_i, 4)) +
                ""
                )

        elif dependent_species == "oh":
            eqn_str_i = ("" +
                "G<sub>OH</sub> = " +
                str(round(slope_i, 4)) +
                " G<sub>OH</sub>+" +
                str(round(intercept_i, 4)) +
                ""
                )

        else:
            eqn_str_i = "TEMP TEMP TEMP TEMP | 190213 | RF"
            raise ValueError('A very specific bad thing happened.')

        annotation_i = dict(
            x=0.,
            y=1.,
            xref="paper",
            yref="paper",
            text=eqn_str_i,
            font=dict(
                color="red",
                family="Droid Sans Mono,Overpass",
                size=9. * (4. / 3.),
                ),
            showarrow=False,
            xanchor="left",
            yshift=-11. * (4. / 3.) * len(self.annotations_list),
            xshift=5.,
            )

        self.annotations_list.append(annotation_i)
        #__|


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

        # print(slope_intercept_dict)

        slope = slope_intercept_dict["slope"]
        intercept = slope_intercept_dict["intercept"]

        def scaling_meth(E_OH):
            """
            """
            #| - scaling_meth
            out = slope * E_OH + intercept

            return(out)
            #__|

        LH_bound = self.x_range[0] - 0.2
        RH_bound = self.x_range[1] + 0.2

        scaling_trace = go.Scatter(
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

        #
        # # Annotation
        # ooh_vs_oh_eqn = ("" +
        #     "G_*OOH = " +
        #     str(round(SC_PLT.scaling_dict["ooh"]["m"], 5)) +
        #     " G_*OH + " +
        #     str(round(SC_PLT.scaling_dict["ooh"]["b"], 5)) +
        #     ""
        #     )
        #
        # o_vs_oh_eqn = ("" +
        #     "G_*O  =  " +
        #     str(round(SC_PLT.scaling_dict["o"]["m"], 5)) +
        #     " G_*OH + " +
        #     str(round(SC_PLT.scaling_dict["o"]["b"], 5)) +
        #     ""
        #     )

        #__|

    # NOTE | This shouldn't be an internal method
    # I don't remember why I wrote the above note
    def get_plotly_layout(self,
        # x_ax_spec="oh",
        title="Scaling Relations",
        showlegend=True,
        layout_dict=None,
        ):
        """Create plotly layout dict.

        Args:
            x_ax_spec:
            title:
            showlegend:
        """
        #| - create_layout
        # if x_ax_spec == ""
        if self.x_ax_species == "oh":
            x_ax_title = "G<sub>ads,*OH</sub> (eV)"
        else:
            x_ax_title = "TEMP"

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

        x_range = self.x_range
        y_range = self.y_range

        layout = {
            "title": None,

            # "titlefont": go.layout.title.Font(size=24),
            # # "titlefont": go.layout.Titlefont(size=24),

            "font": {
                "family": "Arial",  # "Courier New, monospace"
                "color": "black",
                },

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

            #| - Margins ------------------------------------------------------
            "margin": go.layout.Margin(
                b=50.,
                l=50.,
                r=30.,
                t=30.,
                ),
            #__|


            # # Margin
            # "margin": go.layout.Margin(
            #     b=50.,
            #     l=50.,
            #     r=50.,
            #     t=50.,
            #     ),

            "font": dict(
                family='Arial',
                # size=18,
                color='black',
                ),

            #| - Plot Size
            # "width": 37 * 37.795275591,
            # "height": 23 * 37.795275591,
            #
            # # "width": 1.5 * 18.7 * 37.795275591,
            # # "height": 18.7 * 37.795275591,
            # # "width": 1.5 * 18.7 * 37.795275591,
            # # "height": 18.7 * 37.795275591,
            #__|

            #| - Legend
            "showlegend": showlegend,

            "legend": go.layout.Legend(
                x=1.,
                xanchor=None,
                y=1.,
                yanchor="top",
                font=dict(size=10),
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

            # "legend": dict(
            #     font=dict(
            #         size=10,
            #         ),
            #     ),

            #__|

            }

        layout = go.Layout(**layout)

        if layout_dict is not None:
            layout.update(layout_dict)

        # if layout_dict is not None:
        #     from misc_modules.misc_methods import dict_merge
        #     dict_merge(layout_i, layout_dict)
        #     # layout_i = {**layout_i, **layout_dict}

        return(layout)
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
