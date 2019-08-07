#/usr/bin/env python

"""Temporary storage location for methods.

Methods here should be fleshed out and included somewhere in 00_PythonModules
"""

def color_scale_interp(
    input_num,
    max_num,
    min_num,
    color_mesh_size=80,
    hex_mode=True,
    ):
    """
    """
    #| - color_scale_interp
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

    input_norm = ((input_num - min_num)/ (max_num - min_num))

    cs_index = round(input_norm * len(color_scale)) - 1
    if cs_index == -1:
        cs_index = 0
    color_out = color_scale[cs_index]

    if hex_mode:
        def rgb_to_hex(rgb_tuple):
            """
            """
            r = int(rgb_tuple[0])
            g = int(rgb_tuple[1])
            b = int(rgb_tuple[2])

            def clamp(x):
              return max(0, min(x, 255))

            hex_rep = "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))
            return(hex_rep)

        color_out = rgb_to_hex(color_out)


    return(color_out)

    #__|
