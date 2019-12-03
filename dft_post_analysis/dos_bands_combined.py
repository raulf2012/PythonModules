#!/usr/bin/env python

"""Plot pdos and bands side-by-side.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
from plotly import tools
#__|


# pdos_data_out,
# bands_data,
# dos_layout,
# bands_layout,
# plot_title=plot_title,
# subplot_titles=("C-PDOS (N adjacent)", "Bands"),

def plot_pdos_bands(
    pdos_data,
    dos_data,
    bands_data,

    pdos_layout,
    bands_layout,

    e_range=[-15, 10],

    plot_dos=True,
    plot_title="PDOS and Band Diagram",
    subplot_titles=("PDOS", "Bands"),
    ):
    """Plots pdos and bands in adjacent subplots.

    Args:
        pdos_data:
        bands_data:
        pdos_layout:
        bands_layout:
    """
    #| - plot_pdos_bands
    # TEMP_PRINT
    print("KSDKFJDSJIFJSIDJFISD")

    fig = tools.make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True,
        subplot_titles=subplot_titles,
        )

    for trace_i in pdos_data:
        fig.append_trace(trace_i, 1, 1)

    if plot_dos:
        for trace_i in dos_data:
            fig.append_trace(trace_i, 1, 1)

    for trace_i in bands_data:
        fig.append_trace(trace_i, 1, 2)

    fig["layout"].update(
        # height=600,
        # width=600,

        title=plot_title,

        # xaxis1=pdos_layout["xaxis"],
        # xaxis2=bands_layout["xaxis"],

        font={
            "family": "Courier New, monospace",
            "size": 18,
            "color": "black",
            },
        )

    fig["layout"]["yaxis"].update(bands_layout["yaxis"])
    fig["layout"]["yaxis"]["range"] = e_range

    return(fig)
    #__|
