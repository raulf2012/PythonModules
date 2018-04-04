#| - ase_methods - __old__


# def pdos_analysis():
#     """
#     OLD VERSION DON'T USE
#     """

#| - pdos_analysis
# from espresso import espresso
#
# a = read("qn.traj")
#
# calc = espresso(pw=500,             #plane-wave cutoff
#     dw=5000,            #density cutoff
#     xc="BEEF-vdW",      #exchange-correlation functional
#     kpts=(3,3,1),  # ark - k-points for hexagonal symmetry in 2-D materials
#     nbands=-20, #20 extra bands besides the bands needed for valence electrons
#     spinpol = True,     # ark - added spinpolarizatoin
#     sigma=0.1,
#     psppath="/home/vossj/suncat/psp/gbrv1.5pbe",    #pseudopotential path
#     convergence= {
#         "energy": 1.e-5, #convergence parameters
#         "mixing": 0.1,
#         "nmix": 20,
#         "mix": 4,
#         "maxsteps": 500,
#         "diag": "david"
#         },
#     output = {"removesave": False},
#     outdir="calcdir"
#     )    #output directory for Quantum Espresso files
#
# # attach the espresso calculator to the surface
# a.set_calculator(calc)
# a.get_potential_energy()
#
# pdos=calc.calc_pdos(nscf=True,
#                     kpts=(6,6,1),
#                     Emin=-10.0,
#                     Emax=10.0,
#                     tetrahedra=True,
#                     sigma=0.2,
#                     DeltaE=0.01)
#
# pickle.dump(pdos,open("pdos.pkl", "w"))
#__|

#__|
