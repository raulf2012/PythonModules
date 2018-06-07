#!/usr/bin/env python

"""My plot settings and methods.

Author: Raul A. Flores
"""

#| - Import Modules
import plotly

#__|

# "font": {
#     "family": "Courier New, monospace",
#     "size": plot_title_size,
#     "color": "black",
#     },

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
