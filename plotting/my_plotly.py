#!/usr/bin/env python

"""My plot settings and methods.

Author: Raul A. Flores
"""

# | - Import Modules
import os
import copy
import re
import json

import plotly.graph_objs as go
import plotly.io as pl_io
from plotly import io as pyio
# __|

# #########################################################
# | - Nice default plot settings

# | - shared_axis_dict
shared_axis_dict = dict(
    mirror=True,
    showgrid=False,
    showline=True,
    tickcolor='black',
    linecolor='black',
    ticks='outside',  # None, 'outside', 'inside'
    tickfont=dict(
        size=10 * (4 / 3),
        ),
    title=dict(
        font=dict(
            size=12 * (4 / 3),
            ),
        ),
    zeroline=True,
    zerolinecolor='grey',
    zerolinewidth=1.,
    )

shared_axis_dict__bigger = dict(
    mirror=True,
    showgrid=False,
    showline=True,
    tickcolor='black',
    linecolor='black',
    tickfont=dict(
        size=14 * (4 / 3),
        ),
    title=dict(
        font=dict(
            size=16 * (4 / 3),
            ),
        ),
    zeroline=True,
    zerolinecolor='grey',
    zerolinewidth=1.,
    )
# __|

# | - base_plotly_layout
base_plotly_layout = dict(
    paper_bgcolor='rgba(255,255,255,1)',
    plot_bgcolor='rgba(255,255,255,1)',
    font=go.layout.Font(
        color='black',
        family='Arial',
        ),
    yaxis=shared_axis_dict,
    xaxis=shared_axis_dict,
    showlegend=True,
    )
# __|

base_plotly_scatter_props = dict(
    marker={
        'size': 10,
        'line': {
            'width': 0.8,
            'color': 'black',
            },
        }
    )


# __|
# #########################################################


def read_plotly_html(path=None):
    """Read HTML plotly file and return a plotly figure object."""
    #| - read_plotly_html
    fig_path = path

    with open(fig_path) as f:
        html = f.read()
    call_arg_str = re.findall(r'Plotly\.newPlot\((.*)\)', html[-2**16:])[0]
    call_args = json.loads(f'[{call_arg_str}]')
    plotly_json = {'data': call_args[1], 'layout': call_args[2]}
    fig = pl_io.from_json(json.dumps(plotly_json))

    return(fig)
    #__|


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
    tmp_define_both_axis_types=False,
    ):
    """
    Note: range must be set for this to work

    Example usage:

    from plotting.my_plotly import add_duplicate_axes

    shared_axis_data = {
        "tickcolor": "black",
        "ticklen": 3,
        }

    shared_xaxis_data = {
        "dtick": 50,
        **shared_axis_data,
        }

    shared_yaxis_data = {
        "dtick": 1,
        **shared_axis_data,
        }

    shared_meth_props = dict(
        tmp_define_both_axis_types=True,
        )

    add_duplicate_axes(
        fig, axis_type='x',
        axis_data=shared_xaxis_data,
        axis_num_list=[1, ],
        **shared_meth_props)
    add_duplicate_axes(
        fig, axis_type='y',
        axis_data=shared_yaxis_data,
        axis_num_list=[1, ],
        **shared_meth_props)


    """
    # | - add_duplicate_axes

    if axis_type == "x":
        axis_type_other = "y"
    elif axis_type == "y":
        axis_type_other = "x"


    # This is necessary to make sure that the original traces are still visible after adding the new traces
    fig.update_layout(
        # paper_bgcolor="white",
        plot_bgcolor="rgba(255,255,255,0.)",
        )

    axis_info_dict = get_xy_axis_info(fig)[axis_type]
    num_of_axis = axis_info_dict["num_of_axis"]

    # #########################################################################
    if axis_num_list is None:
        axis_num_list = axis_info_dict["axis_num_list"]

    # axis_num_list_new = [i + len(axis_num_list) for i in axis_num_list]
    axis_num_list_new = [i + num_of_axis + 1 for i, j in enumerate(axis_num_list)]

    iterator = enumerate(zip(axis_num_list, axis_num_list_new))
    for i_cnt, (old_index, new_index) in iterator:
        old_Axis = fig.layout[axis_type + "axis" + str(old_index)]

        if old_Axis.range == None:
            print("This doesn't work well if you don't set the range attribute! Do that!")

        new_axis = copy.deepcopy(old_Axis)
        new_axis = new_axis.update(
            showticklabels=False,
            title=dict(
                font=None,
                standoff=None,
                text="",
                ))

        new_axis = new_axis.update(**axis_data)

        axis_key = axis_type + "axis" + str(new_index)
        new_layout = go.Layout({
            axis_key: new_axis,
            })

        fig.update_layout(new_layout)

        if tmp_define_both_axis_types:
            fig.add_scatter(
                **go.Scatter({
                    axis_type + "axis": axis_type + str(new_index),

                    # I added this to fix some issues
                    # It breaks in some applications, look over more closely
                    axis_type_other + "axis": axis_type_other + str(new_index),
                    }).to_plotly_json())

        else:
            fig.add_scatter(
                **go.Scatter({
                    axis_type + "axis": axis_type + str(new_index),
                    }).to_plotly_json())

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

    save_dir=None,
    place_in_out_plot=True,
    plot_name="TEMP_PLOT_NAME",

    write_html=False,
    write_png=False,
    write_pdf=False,
    write_svg=False,
    write_json=True,
    png_scale=6.,  # Not being used now
    # try_orca_write=False,  # COMBAK Not used anymore right?
    # try_kaleido_write=True,
    verbose=False,
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

    #  #  if try_kaleido_write:
    #  import plotly.io as pio
    #  scope = pio.kaleido.scope

    # #####################################################

    if save_dir is None:
        if place_in_out_plot:
            plot_dir = "out_plots"
        else:
            plot_dir = "."
    else:
        if place_in_out_plot:
            plot_dir = os.path.join(save_dir, "out_plots")
        else:
            plot_dir = os.path.join(save_dir, )

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    prepath = os.path.join(plot_dir, plot_name)


    #  # Requires installing the kaleido package to write images
    #  if verbose:
    #      print("Writing pdf with ORCA")
    #      print("prepath:", prepath)

    if write_html:
        fig.write_html(prepath + '.html')

    if write_pdf:
        fig.write_image(prepath + '.pdf')

    if write_svg:
        fig.write_image(prepath + '.svg')

    if write_png:
        fig.write_image(
            prepath + '.png',
            scale=png_scale,
            )

    if write_json:
        fig.write_json(prepath + '.json')


    #| - __old__

    #  # Requires installing the kaleido package to write images
    #  if verbose:
    #      print("Writing pdf with ORCA")
    #      print("prepath:", prepath)

    # -----------------------------------------------------
    # Unused kaleido stuff

    #  import plotly.io as pio
    #  if try_kaleido_write:
    #  scope = pio.kaleido.scope

    #  #  if try_kaleido_write:
    #  scope._shutdown_kaleido()



    # | - Local write to HTML

    # if write_html:
    #     pyio.write_html(
    #         fig,
    #         os.path.join(plot_dir, plot_name + ".html"),
    #         # config=None,
    #         # auto_play=True,
    #         # include_plotlyjs=True,
    #         # include_mathjax=False,
    #         # post_script=None,
    #         # full_html=True,
    #         # animation_opts=None,
    #         # validate=True,
    #         # default_width='100%',
    #         # default_height='100%',
    #         # auto_open=False,
    #         )

    # __|

    # if write_svg:
    #     try:
    #         fig.write_image(prepath + ".svg")
    #     except:
    #         print("Couldn't write svg")

        # try:
        #     # print(prepath)
        #     fig.write_image(prepath + ".png", scale=png_scale)
        # except:
        #     print("Couldn't write png")

    # # Getting the hostname of computer
    # import socket
    # hostname = socket.gethostbyaddr(socket.gethostname())[0]

    # # Requires ORCA installation
    # if (
    #     os.environ.get("USER", "") == "raul-ubuntu-desktop" or
    #     hostname == "raul-ubuntu-vb" or
    #     hostname == "DESKTOP-37GUFJ5" or
    #     hostname == "raul-dell-ubuntu" or
    #     hostname == "raul-dell-latitude" or
    #     try_orca_write
    #     ):

    #__|

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











# | - OLD | add_duplicate_axes

# def plot_layout(
#     # xax_labels =
#     ):
#     """


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
