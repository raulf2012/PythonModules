#!/usr/bin/env python

"""Classes and methods to encapulate colors for plotting.

Author: Raul A. Flores
"""

# | - Import Modules
import seaborn as sns

import colorlover as cl

# __|


# #########################################################
# | - Using seaborn to generate colors
# palette = sns.color_palette(None, 3)
# palette
#
# current_palette = sns.color_palette()
# sns.palplot(current_palette)

def generate_color_palette(
    palette_name="husl",
    bins=10,
    ):
    """

    https://seaborn.pydata.org/tutorial/color_palettes.html

    Args:
        palette_name:
        bins:
            Number of colors to produce
    """
    # | - generate_color_palette
    # tmp = sns.palplot(sns.color_palette(palette_name, bins))
    colors = sns.color_palette(palette_name, bins)

    return(colors)
    # __|

# __|


def color_scale_interp(
    input_num,
    max_num,
    min_num,
    color_mesh_size=80,
    hex_mode=True,
    ):
    """
    """
    # | - color_scale_interp
#     cl.scales["8"]["seq"]["Purples"]

    black_white_cs = [
        'rgb(0,0,0)',
        'rgb(255,255,255)',
        ]

    black_red_cs = [
        'rgb(0,0,0)',
        'rgb(255,0,0)',
        ]

    color_scale_i = black_red_cs

    color_scale_i = cl.scales["8"]["seq"]["Greys"][::-1]
#     color_scale_i = cl.scales["8"]["seq"]["Purples"]
#     color_scale_i = cl.scales['3']['div']['RdYlBu']

    color_scale = cl.interp(
        color_scale_i,
        color_mesh_size,
        )

    color_scale = cl.to_rgb(color_scale)

    # Converting RGB string representatino to tuple
    color_scale = cl.to_numeric(color_scale)

    input_norm = ((input_num - min_num) / (max_num - min_num))

    cs_index = round(input_norm * len(color_scale)) - 1
    if cs_index == -1:
        cs_index = 0

    color_out = color_scale[cs_index]

    if hex_mode:
        color_out = rgb_to_hex(color_out)


    return(color_out)
    # __|
