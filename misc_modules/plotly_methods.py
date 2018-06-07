#!/usr/bin/env python

"""Plotly methods.

Author: Raul A. Flores
"""

#| - Import Modules
import json
# import numpy as np

from plotly.utils import PlotlyJSONEncoder
import plotly.graph_objs as go
#__|

def plotlyfig2json(fig, fpath=None):
    """Serialize a plotly figure object to JSON so it can be persisted to disk.

    Figure's persisted as JSON can be rebuilt using the plotly JSON chart API:

    Code obtained from:
    https://github.com/plotly/plotly.py/issues/579

    http://help.plot.ly/json-chart-schema/

    If `fpath` is provided, JSON is written to file.

    Modified from https://github.com/nteract/nteract/issues/1229
    """
    #| - plotlyfig2json
    redata = json.loads(json.dumps(fig.data, cls=PlotlyJSONEncoder))
    relayout = json.loads(json.dumps(fig.layout, cls=PlotlyJSONEncoder))

    fig_json = json.dumps({'data': redata, 'layout': relayout})

    if fpath:
        with open(fpath, 'w') as f:
            f.write(fig_json)
    else:
        return fig_json
    #__|

def plotlyfromjson(fpath):
    """Render a plotly figure from a json file.

    Code obtained from:
    https://github.com/plotly/plotly.py/issues/579

    Args:
        fpath:
    """
    #| - plotlyfromjson
    with open(fpath, 'r') as f:
        v = json.loads(f.read())

    fig = go.Figure(data=v['data'], layout=v['layout'])

    return(fig)
    # iplot(fig, show_link=False)
    #__|
