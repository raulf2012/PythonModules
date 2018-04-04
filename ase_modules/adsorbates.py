"""Adsorbate related methods."""

#| - IMPORT MODULES
from ase import Atoms
from ase.build import molecule

import numpy as np
#__|

class Adsorbate:
    """Adsorbate atoms object class."""

    #| - Adsorbate ************************************************************

    def __init__(self):
        """Initialize Adsorbate class instance."""
        #| - __init__
        self.tmp = 42
        #__|

    def ooh(self,
        OO_bl=1.359,
        OH_bl=0.993,
        OO_angle=30.0,
        H_up_down="down",
        ):
        """*OOH adsorbate.

        Angles are wrt the z-axis, so 0 is along the z-axis.

        Args:
            OO_bl:
            OO_bl:
            OO_angle:
            H_up_down:
        """
        #| - ooh
        O_1_x = OO_bl * np.sin(np.radians(OO_angle))
        O_1_y = OO_bl * np.cos(np.radians(OO_angle))

        if H_up_down == "up":
            H_angle = 0.0
            H_x = O_1_x + OH_bl * np.sin(np.radians(H_angle))
            H_y = O_1_y + OH_bl * np.cos(np.radians(H_angle))

        elif H_up_down == "down":
            H_angle = 115.0

            H_x = O_1_x + OH_bl * np.sin(np.radians(H_angle))
            H_y = O_1_y + OH_bl * np.cos(np.radians(H_angle))

        ooh_mol = Atoms(["O", "O", "H"],
            positions=[
                (0.0, 0.0, 0.0),
                # (0, 0, 1.359),
                (O_1_x, 0.0, O_1_y),
                # (0.7, 0.0, 2.5),
                (H_x, 0.0, H_y),
                ]
            )

        return(ooh_mol)
        #__|

    def o(self):
        """*O adsorbate."""
        #| - o
        o_mol = Atoms(['O'],
            positions=[
                (0, 0, 0),
                ]
            )

        return(o_mol)
        #__|

    def oh(self,
        OH_bl=0.978,
        ):
        """*OH adsorbate.

        Args:
            OH_bl:
                O-H bond length
        """
        #| - oh
        oh_mol = Atoms(["O", "H"],
            positions=[
                (0, 0, 0),
                (0, 0, OH_bl),
                ]
            )
        return(oh_mol)
        #__|

    def h2o(self,
        H_up_down="up",
        ):
        """H2O adsorbate.

        Args:
            H_up_down:
                Hydrogen atom pointing up or down
        """
        #| - h2o
        h2o_mol = molecule("H2O")

        if H_up_down == "up":
            h2o_mol.rotate(180, "x")
        else:
            pass

        return(h2o_mol)
        #__|

    def get_adsorbate(self, adsorbate, **kwargs):
        """
        Get adsorbate by name.

        Args:
            adsorbate:
            kwargs:
        """
        #| - get_adsorbate
        ads = getattr(self, adsorbate)(**kwargs)

        return(ads)
        #__|
