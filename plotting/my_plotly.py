#!/usr/bin/env python

"""My plot settings and methods.

Author: Raul A. Flores
"""

# | - Import Modules
import os

import copy

# Plotly imports
import plotly

import chart_studio.plotly as py
import plotly.graph_objs as go

from plotly import io as pyio
# __|

def get_xy_axis_info(fig):
    """
    """
    # | - get_xy_axis_info
    xy_axis_list = [i for i in fig.layout._props if "axis" in i]

    x_axis_list = [i for i in xy_axis_list if i[0] == "x"]
    y_axis_list = [i for i in xy_axis_list if i[0] == "y"]

    num_of_xaxis = len(x_axis_list)
    num_of_yaxis = len(y_axis_list)

    # print("x_axis_list:", x_axis_list)
    # print("y_axis_list:", y_axis_list)
    # print("")

    x_axis_num_list = []
    for axis_i in x_axis_list:
        if axis_i == "xaxis" or axis_i == "yaxis":
            x_axis_num_list.append(1)
        else:
            x_axis_num_list.append(int(axis_i[5:]))

    y_axis_num_list = []
    for axis_i in x_axis_list:
        if axis_i == "xaxis" or axis_i == "yaxis":
            y_axis_num_list.append(1)
        else:
            y_axis_num_list.append(int(axis_i[5:]))

    # print("x_axis_num_list:", x_axis_num_list)
    # print("y_axis_num_list:", y_axis_num_list)

    out_dict = dict(
        x=dict(
            num_of_axis=num_of_xaxis,
            axis_list=x_axis_list,
            axis_num_list=x_axis_num_list,
            ),

        y=dict(
            num_of_axis=num_of_yaxis,
            axis_list=y_axis_list,
            axis_num_list=y_axis_num_list,
            ),
        )

    return(out_dict)
    # __|



def add_duplicate_axes(
    fig,
    axis_type="x",  # 'x' or 'y'
    axis_data=dict(),
    axis_num_list=None,
    ):
    """
    """
    # | - add_duplicate_axes

    # This is necessary to make sure that the original traces are still visible after adding the new traces
    fig.update_layout(
        # paper_bgcolor="white",
        plot_bgcolor="rgba(255,255,255,0.)",
        )

    axis_info_dict = get_xy_axis_info(fig)[axis_type]
    num_of_axis = axis_info_dict["num_of_axis"]

    # #########################################################################
    if axis_num_list is None:
        # axis_list = axis_info_dict["axis_list"]
        axis_num_list = axis_info_dict["axis_num_list"]


    # axis_num_list_new = [i + len(axis_num_list) for i in axis_num_list]
    axis_num_list_new = [i + num_of_axis + 1 for i, j in enumerate(axis_num_list)]


    # print("num_of_axis:", num_of_axis)
    # print("axis_num_list_new:", axis_num_list_new)

    # [(i, j) for i, j in enumerate(mylist)]

    iterator = enumerate(zip(axis_num_list, axis_num_list_new))
    for i_cnt, (old_index, new_index) in iterator:
        old_Axis = fig.layout[axis_type + "axis" + str(old_index)]

        new_axis = copy.deepcopy(old_Axis)
        new_axis = new_axis.update(
            # dtick=0.1,
            showticklabels=False,
            title=dict(
                font=None,
                standoff=None,
                text="",
                ),
            )

        new_axis = new_axis.update(**axis_data)

        axis_key = axis_type + "axis" + str(new_index)
        new_layout = go.Layout({
            axis_key: new_axis,
            })

        fig.update_layout(new_layout)

        fig.add_scatter(
            **go.Scatter({
                axis_type + "axis": axis_type + str(new_index)
                }).to_plotly_json())

    # __|


# | - OLD | add_duplicate_axes
# def add_duplicate_axes(
#     fig,
#     axis_type="x",  # 'x' or 'y'
#     axis_data=dict(),
#     ):
#     """
#     """
#     # | - add_duplicate_axes
#
#     # This is necessary to make sure that the original traces are still visible after adding the new traces
#     fig.update_layout(
#         # paper_bgcolor="white",
#         plot_bgcolor="rgba(255,255,255,0.)",
#         )
#
#     # #########################################################################
#     axis_info_dict = get_xy_axis_info(fig)[axis_type]
#
#     num_of_axis = axis_info_dict["num_of_axis"]
#     axis_list = axis_info_dict["axis_list"]
#     axis_num_list = axis_info_dict["axis_num_list"]
#
#     axis_num_list_new = [i + num_of_axis for i in axis_num_list]
#     iterator = enumerate(zip(axis_num_list, axis_num_list_new))
#     for i_cnt, (old_index, new_index) in iterator:
#         old_Axis = fig.layout[axis_type + "axis" + str(old_index)]
#
#         new_axis = copy.deepcopy(old_Axis)
#         new_axis = new_axis.update(
#             # dtick=0.1,
#             showticklabels=False,
#             title=dict(
#                 font=None,
#                 standoff=None,
#                 text="",
#                 ),
#             )
#
#         new_axis = new_axis.update(**axis_data)
#
#         axis_key = axis_type + "axis" + str(new_index)
#         new_layout = go.Layout({
#             axis_key: new_axis,
#             })
#
#         fig.update_layout(new_layout)
#
#         fig.add_scatter(
#             **go.Scatter({
#                 axis_type + "axis": axis_type + str(new_index)
#                 }).to_plotly_json())
#     # __|
# __|


def add_minor_ticks(
    fig,
    axis="x",  # 'x', 'y', or 'both'
    ticks_props_new_x=None,
    ticks_props_new_y=None,
    ):
    """
    """
    # | - add_minor_ticks
    dummy_trace = go.Scatter(
        x=[0],
        y=[0],
        # xaxis="x2",
        xaxis="x",
        yaxis="y",
        opacity=0.,
        name="TEMP|8asdf",
        )

    xaxis_orig = copy.deepcopy(fig.layout.xaxis)
    yaxis_orig = copy.deepcopy(fig.layout.yaxis)


    if axis == "x" or axis == "both":
        for trace in fig.data:
            if trace.xaxis == None or trace.xaxis == "x":
                xaxis_new = "x2"
            else:
                tmp = "x2"
                xaxis_old = trace.xaxis
                tmp = int(xaxis_old[1:])
                xaxis_new = "x" + str(tmp)

            fig.add_scatter(
                **trace.update(dict1=dict(xaxis=xaxis_new), overwrite=True).to_plotly_json(),
                )





    if axis == "y" or axis == "both":
        for trace in fig.data:
            if trace.yaxis == None or trace.yaxis == "y":
                yaxis_new = "y2"
            else:
                tmp = "y2"
                yaxis_old = trace.yaxis
                tmp = int(axis_old[1:])
                yaxis_new = "y" + str(tmp)

            fig.add_scatter(
                **trace.update(dict1=dict(yaxis=yaxis_new), overwrite=True).to_plotly_json(),
                )




    # #########################################################################
    # global_axis_props = go.layout.XAxis(

    global_axis_props = dict(
        showticklabels=False,
        # title=go.layout.xaxis.Title(
        title=dict(
            font=None,
            standoff=None,
            # text="d",
            text="",
            ),
        )

    new_xaxis = xaxis_orig.update(
        # global_axis_props.to_plotly_json(),
        global_axis_props,
        overwrite=True)
    new_xaxis = new_xaxis.update(ticks_props_new_x)

    new_yaxis = yaxis_orig.update(
        global_axis_props,
        overwrite=True)
    # new_yaxis = new_yaxis.update(ticks_props_new_y.to_plotly_json())
    new_yaxis = new_yaxis.update(ticks_props_new_y)

    tmp = fig.update_layout(
        xaxis2=new_xaxis,
        yaxis2=new_yaxis,
        )



    fig.add_scatter(**dummy_trace.to_plotly_json())

    return(fig)

    # __|


def my_plotly_plot(
    figure=None,
    plot_name="TEMP_PLOT_NAME",
    write_html=False,
    write_png=False,
    png_scale=6.,
    write_pdf=False,
    write_svg=False,

    try_orca_write=False,
    ):
    """
    Returns: Plotly figure object


    TODO Add functionality to display image using Ipython.display instead of displaying the full plotly html (interactive), this will save space

    from IPython.display import Image
    Image("out_plot/" + plot_name_i + ".png")

    Args:
    ---------------------------------------------------------------------------
    layout:
      plotly layout
    layout_override:
      Dictionary to override layout
    plot_name:
      Plot name (used both for plot upload and local save)
    data:
      plotly data object

    """
    # | - my_plotly_plot
    assert figure is not None, "Must pass a plot.ly figure object"
    fig = figure


    # #########################################################################
    plot_dir = "out_plot"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    # | - Local write to HTML
    if write_html:
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
    # __|

    # | - Write pdf and svg (if ORCA is installed and working)
    # Getting the hostname of computer
    import socket
    hostname = socket.gethostbyaddr(socket.gethostname())[0]

    # Requires ORCA installation
    if (
        os.environ.get("USER", "") == "raul-ubuntu-desktop" or
        hostname == "raul-ubuntu-vb" or
        hostname == "DESKTOP-37GUFJ5" or
        hostname == "raul-dell-ubuntu" or
        hostname == "raul-dell-latitude" or
        try_orca_write
        ):
        print("Writing pdf with ORCA")

        prepath = os.path.join(plot_dir, plot_name)
        print("prepath:", prepath)

        if write_pdf:
            try:
                fig.write_image(prepath + ".pdf")
            except:
                print("Couldn't write pdf")
        if write_svg:
            try:
                fig.write_image(prepath + ".svg")
            except:
                print("Couldn't write svg")
        if write_png:
            try:
                fig.write_image(prepath + ".png", scale=png_scale)
            except:
                print("Couldn't write png")

    # __|


    # return(fig)
    # __|




def reapply_colors(data):
    """Redefines the line colors of a plotly data series.

    Groups by legend grouping, fix this it's not general enough

    Args:
        plotly data series (list of graph objects to be plotted)
    """
    # | - reapply_colors
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
    # __|

def plot_layout(
    # xax_labels =
    ):
    """
    """
    # | - plot_layout

    # | - Plot Settings
    plot_title_size = 18
    tick_lab_size = 16
    axes_lab_size = 18
    legend_size = 18
    # __|

    # | - Plot Layout
    xax_labels = ["O2", "OOH", "O", "OH", "H2O"]
    layout = {

        "title": "FED of ORR Mechanism For Iron-Supported-Graphene",

        "font": {
            "family": "Courier New, monospace",
            "size": plot_title_size,
            "color": "black",
            },

        # | - Axes --------------------------------------------------------------
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
        # __| -------------------------------------------------------------------

        # | - Legend ------------------------------------------------------------
        "legend": {
            "traceorder": "normal",
            "font": dict(size=legend_size)
            },
        # __| -------------------------------------------------------------------

        # | - Plot Size
        "width": 200 * 4.,
        "height": 200 * 3.,
        # __|

        }
    # __|

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

    # __|
