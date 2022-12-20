#!/usr/bin/env python

"""Classes and methods to encapulate colors for plotting.

Author: Raul A. Flores
"""

# | - Import Modules
import matplotlib.pyplot as plt
import random

import matplotlib.colors
# __|


# #########################################################
# | - My Sets of Contrasting Color Schemes

color_list = [
    "rgb(113,166,190)",
    "rgb(145,209,79)",
    "rgb(124,78,196)",
    "rgb(203,169,87)",
    "rgb(200,72,150)",
    "rgb(130,204,159)",
    "rgb(190,82,63)",
    "rgb(80,51,82)",
    "rgb(81,92,54)",
    "rgb(192,151,188)",
    ]

color_palette_1 = [
    "rgb(198,158,61)",
    "rgb(238,88,121)",
    "rgb(236,186,51)",
    "rgb(237,101,67)",
    "rgb(222,134,54)",
    ]

color_palette_2 = [
    "rgb(181,149,213)",
    "rgb(103,76,204)",
    "rgb(97,166,201)",
    "rgb(185,76,198)",
    "rgb(93,80,139)",
    ]

# __|

def rgb_to_hex(rgb_tuple):
    """
    """
    # | - rgb_to_hex
    r = int(rgb_tuple[0])
    g = int(rgb_tuple[1])
    b = int(rgb_tuple[2])

    def clamp(x):
        return(max(0, min(x, 255)))

    hex_rep = "#{0:02x}{1:02x}{2:02x}".format(
        clamp(r),
        clamp(g),
        clamp(b),
        )

    return(hex_rep)
    # __|

def get_random_color():
    """Generate random color in hex format."""
    # | - get_random_color
    number_of_colors = 1

    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                 for i in range(number_of_colors)]

    return(color[0])
    # __|

def rgb_to_hex(rgb_tuple):
    """
    """
    # | - rgb_to_hex
    r = int(rgb_tuple[0])
    g = int(rgb_tuple[1])
    b = int(rgb_tuple[2])

    def clamp(x):
        return(max(0, min(x, 255)))

    hex_rep = "#{0:02x}{1:02x}{2:02x}".format(
        clamp(r),
        clamp(g),
        clamp(b),
        )

    return(hex_rep)
    # __|

def hex_to_rgb(hex_code):
    """
    """
    # | - hex_to_rgb
    return(
        matplotlib.colors.to_rgb(hex_code))
    #  __|

def darken_color(color, darken_amount=0.2):
    """
    """
    # | - darken_color
    new_color_tuple = []
    for i in tuple([i - darken_amount for i in color]):
        if i < 0:
            new_color_tuple.append(0)
        else:
            new_color_tuple.append(i)
    color_darkened = rgb_to_hex(new_color_tuple)

    return(color_darkened)
    # __|

def lighten_color(color, lighten_amount=0.2):
    """
    """
    # | - lighten_color
    new_color_tuple = []
    for i in tuple([i + lighten_amount for i in color]):
        if i > 1:
            new_color_tuple.append(1)
        else:
            new_color_tuple.append(i)
    color_lightened = rgb_to_hex(new_color_tuple)

    return(color_lightened)
    # __|
