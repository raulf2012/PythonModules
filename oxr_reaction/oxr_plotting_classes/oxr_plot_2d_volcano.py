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


# ██    ██  ██████  ██       ██████         ██████  ██       ██████  ████████
# ██    ██ ██    ██ ██      ██              ██   ██ ██      ██    ██    ██
# ██    ██ ██    ██ ██      ██              ██████  ██      ██    ██    ██
#  ██  ██  ██    ██ ██      ██              ██      ██      ██    ██    ██
#   ████    ██████  ███████  ██████ ███████ ██      ███████  ██████     ██

class Volcano_Plot_2D():
    """Class to plot OER/ORR volcano plots.

    Development Notes:
        TEMP
    """

    #| - Volcano_Plot *********************************************************

    def __init__(self,
        ORR_Free_E_Plot,
        plot_range=None,
        ):
        """
        """
        #| - __init__
        self.ORR_Free_E_Plot = ORR_Free_E_Plot

        if plot_range is None:
            self.plot_range = {
                "x": [+0.5, +2.5],
                "y": [-0.5, +2.0],
                }
        else:
            self.plot_range = plot_range

        # #####################################################################


        self.contour_trace = self.__create_contour_trace__()
        self.data_point_traces = self.__create_data_point_traces__()


        data = []
        data += self.data_point_traces
        data.insert(0, self.contour_trace)
        self.traces = data
        #__|

    def __create_contour_trace__(self):
        """
        """
        #| - __create_contour_trace__
        x_range_bounds = self.plot_range["x"]
        y_range_bounds = self.plot_range["y"]


        x_range = np.arange(
            x_range_bounds[0] - 1.,
            x_range_bounds[1] + 1.,
            0.01)

        y_range = np.arange(
            y_range_bounds[0] - 1.,
            y_range_bounds[1] + 1.,
            0.01)

        # x_range = np.arange(0.9, 2, 0.01)
        # y_range = np.arange(-0.5, 2., 0.01)

        X, Y = np.meshgrid(x_range, y_range)

        for ind_i, x_i in enumerate(x_range):
            for ind_j, oh_j in enumerate(y_range):
                X[ind_j][ind_i] = self.overpotential3(x_i, oh_j)

        series_i = go.Contour(
            z=X,
            x=x_range,
            y=y_range,

            zmin=0.2,
            zmax=1.2,
        #     zmax=1.2,

        #     ncontours=180,
            ncontours=20,

            #| - Color Scale
            colorscale='Jet',
            reversescale=True,
            autocontour=True,

            # showscale=False,
            #__|

            #| - Colorbar
            colorbar=go.contour.ColorBar(
                x=None,
                xanchor=None,
                xpad=None,
                y=None,
                yanchor=None,
                ),
            #__|

            #| - Line
            line=go.contour.Line(
                color="white",
                # dash="dot",
                smoothing=1.,
                width=0.5,
                )
            )
            #__|

        return(series_i)
        #__|



    def __create_data_point_traces__(self):
        """
        """
        #| - __create_data_point_traces__
        data_list = []
        for sys_i in self.ORR_Free_E_Plot.series_list:
            trace_i = self.__create_scatter_trace_i__(sys_i)
            data_list.append(trace_i)

        return(data_list)
        #__|

    def __create_scatter_trace_i__(self,
        sys_i,
        smart_format_i=None,
        ):
        """
        """
        #| - __create_trace_i__
        trace_i = go.Scatter(
            x=[sys_i.energy_states_dict["o"] - sys_i.energy_states_dict["oh"]],
            y=[sys_i.energy_states_dict["oh"]],
            mode='markers',

            # name=sys_i.series_name,
            name=sys_i.series_name,

            # <br>

            marker=dict(
                size=20,
                symbol=sys_i.format_dict["symbol_i"],
                color=sys_i.format_dict["color_2"],
                line=dict(
                    width=2,
                    color=sys_i.format_dict["color_1"],
                    ),
                ),
            )

        return(trace_i)
        #__|










    # #########################################################################
    # #########################################################################

    def ooh_oh_scaling(self, doh):
        """ooh_oh_scaling equation."""
        #| - ooh_oh_scaling
        #like ambars
        #dooh=0.5*doh  + 3.0		 #O
        #normal one
        # dooh = doh + 3.2

        # m = 0.8926
        # b = 3.174

        m = 0.9104
        b = 3.144

        # 0.9104 x + 3.144
        # 1.302 x + 1.338

        dooh = m * doh + b  # TEMP

        return(dooh)
        #__|

    def overpotential3(self, x, doh):
        """Calculate overpotential (version 3).

        Args:
            x:
            doh:
        """
        #| - overpotential3
        dooh = self.ooh_oh_scaling(doh)
        dg14 = [doh, x, dooh - (x + doh), -dooh + 4.92]
        m = max(dg14)
        return(m - 1.23)

        #return doh*do
        #__|

    # #########################################################################
    # #########################################################################


    def get_plotly_layout(self,
        showlegend=False,
        layout_dict=None,

        # height=9. * 37.795275591,
        # width=9. * 37.795275591,
        ):
        """
        """
        #| - get_plotly_layout
        y_range = self.plot_range["y"]
        x_range = self.plot_range["x"]

        tick_lab_size = 8 * (4. / 3.)
        axes_lab_size = 9 * (4. / 3.)

        layout = {

            "title": None,

            #| - Font Settings
            "font": {
                "family": "Arial",  # "Courier New, monospace"
                "color": "black",
                },
            #__|

            #| - Axes ---------------------------------------------------------

            #| - yaxis
            "yaxis": {
                "title": "ΔG<sub>OH</sub> (eV)",

                # "title": "$\\Delta G (ev)$",
                "range": y_range,

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
                "ticks": 'outside',

                # "tick0": 0,
                # "dtick": 0.1,

                "tickcolor": 'black',
                "ticklen": 2,
                "tickwidth": 1,
                },
            #__|

            #| - xaxis
            "xaxis": {
                "title": "ΔG<sub>O</sub> - ΔG<sub>OH</sub> (eV)",
                "range": x_range,

                "zeroline": False,
                "showline": True,
                "mirror": True,
                # "linecolor": 'red',
                "linecolor": 'black',
                "showgrid": False,
                "titlefont": dict(size=axes_lab_size),
                "showticklabels": True,
                "ticks": 'outside',

                # "tick0": 0,
                # "dtick": 0.2,

                "ticklen": 2,
                "tickwidth": 1,
                "tickcolor": 'black',
                "tickfont": dict(
                    size=tick_lab_size,
                    ),
                },
            #__|

            #__|

            #| - Margins ------------------------------------------------------
            "margin": go.layout.Margin(
                b=50.,
                l=50.,
                r=30.,
                t=30.,
                ),
            #__|

            #| - Legend ---------------------------------------------------
            "legend": go.layout.Legend(
                x=1.2,
                xanchor=None,
                y=0.98,
                yanchor="top",
                font=dict(size=18),
                bgcolor="rgba(200,255,255,0.85)",

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

            "showlegend": True,
            # "showlegend": showlegend,
            #__|

            #| - Plot Size
            "width": 37 * 37.795275591,
            "height": 23 * 37.795275591,
            #__|

            "paper_bgcolor": 'rgba(255,255,255,1.)',
            # "plot_bgcolor": 'rgba(3,3,3,0.3)',
            # "plot_bgcolor": 'rgba(0,0,0,0)',

            }

        #| - Plot Size Settings
        # # bottom_margin_size = 2.5 * 9. * 37.795275591
        # plot_size_settings = {
        #
        #     # "width": 37 * 37.795275591,
        #     # "height": 23 * 37.795275591,
        #
        #     # "width": 25 * 37.795275591,
        #     # "height": 17 * 37.795275591,
        #     }
        # layout = {**layout, **plot_size_settings}
        #__|

        layout = go.Layout(**layout)

        if layout_dict is not None:
            layout.update(layout_dict)

            # from misc_modules.misc_methods import dict_merge
            # dict_merge(layout, layout_dict)


        return(layout)  # COMBAK

        #__|

    #__| **********************************************************************