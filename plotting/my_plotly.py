#!/usr/bin/env python

"""My plot settings and methods.

Author: Raul A. Flores
"""

#| - Import Modules
import plotly

import os
# import plotly.plotly as py
import chart_studio.plotly as py
import plotly.graph_objs as go

from plotly import io as pyio
#__|


#| - Plotly
def reapply_colors(data):
    """Redefines the line colors of a plotly data series.

    Groups by legend grouping, fix this it's not general enough

    Args:
        plotly data series (list of graph objects to be plotted)
    """
    #| - reapply_colors
    from colors.colors import generate_color_palette

    dat_lst_master = data

    groups_list = []
    for series_i in dat_lst_master:
        groups_list.append(series_i.legendgroup)

    groups_list = list(set(groups_list))
    # print(groups_list)

    num_series = len(groups_list)
    colors = generate_color_palette(bins=num_series)

    new_list = []
    for color in colors:
        color_new = tuple([int(255 * x) for x in color])
        color_new = "rgb" + str(color_new)
        new_list.append(color_new.replace(" ", ""))
    colors = new_list

    colors_dict = dict(zip(groups_list, colors))

    for series_i in dat_lst_master:
        tmp = colors_dict[series_i.legendgroup]

        series_i.marker["color"] = tmp
        series_i.line["color"] = tmp

    return(dat_lst_master)
    #__|

def plot_layout(
    # xax_labels =
    ):
    """
    """
    #| - plot_layout

    #| - Plot Settings
    plot_title_size = 18
    tick_lab_size = 16
    axes_lab_size = 18
    legend_size = 18
    #__|

    #| - Plot Layout
    xax_labels = ["O2", "OOH", "O", "OH", "H2O"]
    layout = {

        "title": "FED of ORR Mechanism For Iron-Supported-Graphene",

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
        "width": 200 * 4.,
        "height": 200 * 3.,
        #__|

        }
    #__|

    fig = Figure(data=dat_lst, layout=layout)
    # plotly.plotly.image.save_as(fig, filename="pl_hab_opda_raman.png")

    plotly.offline.plot(
        {
            "data": dat_lst,
            "layout": layout,
            },
        filename="plots/pl_fed_supp_graph_02.html"
        )

    # tmp = plotly.plotly.image.plot(data, filename="pl_fed_180314.png")

    return(layout)

    #__|

#__|



def my_plotly_plot(
    figure=None,
    layout=None,
    layout_override=None,
    plot_name=None,
    save_dir=None,
    data=None,
    upload_plot=True,
    ):
    """
    TODO:
      Remove layout override functionality, this should be done before calling
      the method

    Returns: Plotly figure object

    Args:
    ---------------------------------------------------------------------------
    layout:
      plotly layout
    layout_override:
      Dictionary to override layout
    plot_name:
      Plot name (used both for plot upload and local save)
    save_dir:
      Plot.ly folder to save figure into (Not used for local save)
    data:
      plotly data object
    upload_plot:
      Upload plot to plotly servers

    """
    #| - plot_proto
    if layout is None:
        layout = go.Layout()

    if figure is not None:
        fig = figure
    else:
        fig = go.Figure(data=data, layout=layout)


    fig.layout.update(layout_override)


    #| - Upload to plot.ly website
    # #########################################################################
    if upload_plot:
        plotly_filename = os.path.join(
            save_dir,
            # "02_oer_analysis",
            # "oer_2d_volcano_plot",
            plot_name)
        tmp = py.iplot(fig, filename=plotly_filename)
        print(plotly_filename)
    #__|


    # #########################################################################
    plot_dir = "out_plot"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    #| - Local write to HTML
    pyio.write_html(
        fig,
        os.path.join(plot_dir, plot_name + ".html"),
        # config=None,
        # auto_play=True,
        # include_plotlyjs=True,
        # include_mathjax=False,
        # post_script=None,
        # full_html=True,
        # animation_opts=None,
        # validate=True,
        # default_width='100%',
        # default_height='100%',
        # auto_open=False,
        )
    #__|


    #| - Write pdf and svg (if ORCA is installed and working)
    import socket
    hostname = socket.gethostbyaddr(socket.gethostname())[0]

    # Requires ORCA installation
    if os.environ["USER"] == "raul-ubuntu-desktop" or hostname == "raul-ubuntu-vb":
        print("Writing pdf with ORCA")

        # This seems to be the preferred syntax now
        fig.write_image(
            os.path.join(plot_dir, plot_name + ".pdf")
            # "out_plot/test_fig.pdf"
            )

        # This seems to be the preferred syntax now
        fig.write_image(
            os.path.join(plot_dir, plot_name + ".svg")
            # "out_plot/test_fig.pdf"
            )
    #__|


    return(fig)
    #__|
