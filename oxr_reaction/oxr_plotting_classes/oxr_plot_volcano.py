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
# from orr_reaction.orr_series import ORR_Free_E_Series

from oxr_reaction.adsorbate_scaling import lim_U_i
# from orr_reaction.adsorbate_scaling import lim_U_i
#__|


# ██    ██  ██████  ██       ██████         ██████  ██       ██████  ████████
# ██    ██ ██    ██ ██      ██              ██   ██ ██      ██    ██    ██
# ██    ██ ██    ██ ██      ██              ██████  ██      ██    ██    ██
#  ██  ██  ██    ██ ██      ██              ██      ██      ██    ██    ██
#   ████    ██████  ███████  ██████ ███████ ██      ███████  ██████     ██

class Volcano_Plot():
    """Class to plot OER/ORR volcano plots.

    Development Notes:
        TEMP
    """

    #| - Volcano_Plot *********************************************************

    def __init__(self,
        ORR_Free_E_Plot,
        x_ax_species="o-oh",  # 'o-oh' or 'oh'
        smart_format_dict=None,
        plot_range=None,
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
            plot_range:
                Ex.)
                plot_range={
                    "y": [1., 5.],
                    "x": [-2., 4.],
                    }

        """
        #| - __init__
        self.ORR_Free_E_Plot = ORR_Free_E_Plot
        self.x_ax_species = x_ax_species
        self.plot_range = plot_range
        self.smart_format_dict = smart_format_dict

        self.data_points = []
        self.data_lines = []

        self.marker_color_key = marker_color_key
        self.marker_border_color_key = marker_border_color_key
        self.marker_shape_key = marker_shape_key

        #__|

    # NOTE | Rename this create_volcano_plot
    def create_volcano_relations_plot(self,
        show_data_labels=False,
        # smart_format_dict=None,
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
                [{"bulk_system": "IrO3"}, {self.marker_color_key: "green"}],
                [{"bulk_system": "IrO2"}, {self.marker_color_key: "yellow"}],

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

            name_i = series_i.series_name

            if series_i.color is not None:
                smart_format_i[self.marker_color_key] = series_i.color


            format_i = smart_format_i

            if series_i.format_dict:
                format_i = series_i.format_dict

            trace_i = self.__create_trace_i__(
                x_ax_energy,
                lim_pot_i,
                # smart_format_i,
                format_i,
                name_i,
                group=series_i.group,
                show_data_labels=show_data_labels,
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
        trace_priority="top",  # 'top' or 'bottom'
        legs_to_plot=[
            "o2_to_ooh",
            "ooh_to_o",
            "o_to_oh",
            "oh_to_h2o",
            ],
        line_color="black",
        ):
        """Create volcano data traces.

        Args:
            gas_molec_dict:
            scaling_dict:
            plot_all_legs:
            plot_min_max_legs:
            trace_priority:
                if 'top', the volcano lines will be placed on the top of the
                plot, if 'bottom' then the data points will by on top of the
                volcano
        """
        #| - create_volcano_lines
        out_data = []
        x_range = self.plot_range["x"]

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

        #| - Create Volcano Legs (LOOP)
        x_axis = np.linspace(x_range[0], x_range[1], num=500)
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
        #__|

        if plot_all_legs:

            #| - plot_all_legs
            # hoverinfo_type = "none"
            hoverinfo_type = "name"

            trace_o2_to_ooh = go.Scatter(
                x=x_axis,
                y=energy_dict["o2_to_ooh"],
                name="O2->*OOH",
                hoverinfo=hoverinfo_type,
                line=dict(
                    color="#e7b8bc",
                    width=2,
                    dash="solid",
                    )
                )

            trace_ooh_to_o = go.Scatter(
                x=x_axis,
                y=energy_dict["ooh_to_o"],
                name="*OOH->*O",
                hoverinfo=hoverinfo_type,
                line=dict(
                    color="#afd7c3",
                    width=2,
                    dash="solid",
                    )
                )

            trace_o_to_oh = go.Scatter(
                x=x_axis,
                y=energy_dict["o_to_oh"],
                name="*O->*OH",
                hoverinfo=hoverinfo_type,
                line=dict(
                    color="#b5c4e2",
                    width=2,
                    dash="solid",
                    )
                )

            trace_oh_to_h2o = go.Scatter(
                x=x_axis,
                y=energy_dict["oh_to_h2o"],
                name="*OH->H2O",
                hoverinfo=hoverinfo_type,
                line=dict(
                    color="#dbcdab",
                    width=2,
                    dash="solid",
                    )
                )

            if trace_priority == "top":
                out_data.append(trace_o2_to_ooh)
                out_data.append(trace_ooh_to_o)
                out_data.append(trace_o_to_oh)
                out_data.append(trace_oh_to_h2o)

            elif trace_priority == "bottom":
                out_data.insert(0, trace_o2_to_ooh)
                out_data.insert(0, trace_ooh_to_o)
                out_data.insert(0, trace_o_to_oh)
                out_data.insert(0, trace_oh_to_h2o)
            #__|

        #__|

        #| - Minimum Energy Legs
        energy_lists= []
        for leg_i in legs_to_plot:
            energy_lists.append(energy_dict[leg_i])

        min_max_e_list = []
        for legs in zip(*energy_lists):
            if self.ORR_Free_E_Plot.rxn_type == "ORR":
                energy_i = min(*legs)

            elif self.ORR_Free_E_Plot.rxn_type == "OER":
                energy_i = max(*legs)

            min_max_e_list.append(energy_i)

        trace_volcano = go.Scatter(
            x=x_axis,
            y=min_max_e_list,
            name="activity volcano",
            hoverinfo="skip",
            line=dict(
                color=line_color,
                width=2,
                # dash="dash",
                dash="5px,2px,5px,2px",
                )
            )

        if plot_min_max_legs:
            if trace_priority == "top":
                out_data.append(trace_volcano)

            elif trace_priority == "bottom":
                out_data.insert(0, trace_volcano)
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
        show_data_labels=False,
        ):
        """
        """
        #| - __create_trace_i__

        if show_data_labels is True:
            mode_i = "markers+text"
        elif show_data_labels is False:
            mode_i = "markers"
        else:
            print("TEMP TEMP TEMP | __create_trace_i__")

        # print(mode_i)

        trace_i = go.Scatter(
            x=[x_energy],
            y=[y_energy],

            mode=mode_i,
            # mode="markers+text",
            # mode="markers",

            name=name_i,
            text=[name_i],
            # text=["TEMP"],

            legendgroup=group,

            hoverinfo="text",

            # hoverinfo=None,
            # hoverinfosrc=None,
            # hoverlabel=None,
            # hoveron=None,
            # hovertext=None,
            # hovertextsrc=None,

            # textposition='top right',
            textposition='middle left',
            textfont={
                # "family": "Courier New, monospace",
                # "family": font_family,
                "size": 10,
                "color": "black",
                },

            marker=dict(
                size=smart_format_i.get("marker_size", 9),
                color=smart_format_i.get(self.marker_color_key, "red"),
                symbol=smart_format_i.get(
                    self.marker_shape_key, "circle"),
                line=dict(
                    width=smart_format_i.get("marker_border_width", 1.),
                    color=smart_format_i.get(
                        self.marker_border_color_key, "black"),
                    ),
                ),
            )

        return(trace_i)
        #__|

    def get_plotly_layout(self,
        showlegend=False,
        width=9. * 37.795275591,
        height=9. * 37.795275591,
        layout_dict=None,
        ):
        """
        """
        #| - get_plotly_layout

        #| - Properties
        # plot_title="FED"
        plot_title = None
        # plot_title_size = 18
        # tick_lab_size = 9 * (4. / 3.)
        tick_lab_size = 8 * (4. / 3.)
        axes_lab_size = 9 * (4. / 3.)
        legend_size = 18
        # font_family="Computer Modern"  # "Courier New, monospace"
        font_family = "Arial"  # "Courier New, monospace"
        #__|

        # self.x_ax_spec

        if self.x_ax_species == "oh":
            # xaxis_title = "dG_*OH (eV)"
            xaxis_title = "dG<sub>OH</sub> (eV)"
        elif self.x_ax_species == "o-oh":
            # xaxis_title = "dG_*OH - dG_*O (eV)"
            xaxis_title = "dG<sub>O</sub> - dG<sub>OH</sub> (eV)"

        # layout["margin"] = go.layout.Margin(
        #     b=50.,
        #     l=50.,
        #     r=50.,
        #     t=50.,
        # #     pad=20.,
        #     )

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
                # "linecolor": 'red',
                "linecolor": 'black',
                "showgrid": False,

                "titlefont": dict(size=axes_lab_size),

                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                "ticks": 'inside',
                "tick0": 0,
                "tickcolor": 'black',
                "dtick": 0.1,
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
                # "linecolor": 'red',
                "linecolor": 'black',
                "showgrid": False,
                "titlefont": dict(size=axes_lab_size),
                "showticklabels": True,
                "ticks": 'inside',
                "tick0": 0,
                "dtick": 0.2,
                "ticklen": 2,
                "tickwidth": 1,
                "tickcolor": 'black',
                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                },
            #__|

            #__|

            "margin": go.layout.Margin(
                b=50.,
                l=50.,
                r=50.,
                t=50.,
                ),

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

        #| - Applying Layout override dict
        if layout_dict is not None:
            from misc_modules.misc_methods import dict_merge
            dict_merge(layout, layout_dict)

            # layout_i = {**layout_i, **layout_dict}

        #__|

        return(layout)

        #__|


    #__| **********************************************************************
