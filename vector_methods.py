"""Methods to compute properties of vectors.

Author: Raul A. Flores

Development Notes:
    TEMP
"""


# | - Import Modules
import os

import numpy as np
# __|



def radians_to_degrees(angle):
    """
    """
    # | - radians_to_degrees
    angle_degrees = angle * 180 / np.pi

    return(angle_degrees)
    # __|

def degrees_to_radians(angle):
    """
    """
    # | - radians_to_degrees
    angle_radians = angle * (180 / np.pi) ** -1

    return(angle_radians)
    # __|

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    # | - unit_vector
    return vector / np.linalg.norm(vector)
    # __|

def angle_between_vectors(v1, v2, return_radians=True):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    # | - angle_between_vectors
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if not return_radians:
        angle = radians_to_degrees(angle)

    return(angle)
    # __|

def get_vector_magnitude(v1):
    # | - get_vector_magnitude
    magnitude = np.linalg.norm(v1)
    return(magnitude)
    # __|

def get_rotation_matrix(theta, axis):
    """
    https://en.wikipedia.org/wiki/Rotation_matrix
    """
    # | - get_get_rotation_matrix

    rotation_matrix_x = np.array([
        [1., 0.,             0.,            ],
        [0., np.cos(theta), -np.sin(theta), ],
        [0., np.sin(theta),  np.cos(theta), ],
        ])

    rotation_matrix_y = np.array([
        [ np.cos(theta), 0., np.sin(theta), ],
        [ 0.,            1., 0.,            ],
        [-np.sin(theta), 0., np.cos(theta), ],
        ])

    rotation_matrix_z = np.array([
        [np.cos(theta), -np.sin(theta), 0., ],
        [np.sin(theta),  np.cos(theta), 0., ],
        [0.,             0.,            1., ],
        ])

    if axis == "x":
        rotation_matrix = rotation_matrix_x
    if axis == "y":
        rotation_matrix = rotation_matrix_y
    if axis == "z":
        rotation_matrix = rotation_matrix_z

    return(rotation_matrix)
    # __|
