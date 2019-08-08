#| - IMPORT MODULES
import numpy as np
#__|

def unit_vector(vector):
    """ Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """
    #| - unit_vector
    return vector / np.linalg.norm(vector)
    #__|

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """
    #| - angle_between
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    #__|

def smooth_data_series(y_data, box_pts):
    """Smooth data series by convolution
    """
    #| - smooth_data_series
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y_data, box, mode="same")
    return(y_smooth)
    #__|

def make_filter_list(len_data, percent_keep):
    """
    """
    #| - filter_list

#     len_data = len(pdos_data[0])

    filter_list = np.random.choice(
        [True, False],
        size=len_data,
        p=[percent_keep, 1. - percent_keep]
        )

    return(filter_list)
    #__|
