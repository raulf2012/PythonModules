#!/usr/bin/env python

"""Classes and methods to encapulate colors for plotting.

Author: Raul A. Flores
"""

#| - Import Modules
import seaborn as sns
#__|

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
    #| - generate_color_palette
    # tmp = sns.palplot(sns.color_palette(palette_name, bins))
    colors = sns.color_palette(palette_name, bins)

    return(colors)
    #__|
