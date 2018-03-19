
#| - __old__

# def process_orr_adsorbates_fed(free_e_dict):
#     """
#     Calculates the overpotential for the ORR reaction given the Free Energies
#     of the species. Also outputs the Free energy landscape at V=0 and V=V_eq,
#     as well as the x-coordinate data (just integer pairs defining the steps).
#
#     Output: {'overpotential': ###, 'ydata_nobias': [##,##,...], 'ydata_eq':
#     [##,##,...], 'xdata': [0,1,1,2,2,...]}
#
#     Args
#         free_e_lst: Dictionary containing the free energies of intermediates
#         format the ORR reaction. An example is shown below.
#         Ex.
#         {'Ni': {'bulk': 0, 'ooh': 4.87, 'o': 4.17, 'oh': 2.11}}
#     """
#     #| - process_orr_adsorbates_fed
#
#     #| - Appending Species Free Energies to Output Dict                   1111
#     metal = free_e_dict.keys()[0]
#     species_fe = free_e_dict[metal]
#     #__|                                                                  1111
#
#     fe_all = free_e_dict
#
#     steps = range(5)
#
#     #| - Constructing Reaction Coordinate - AUTOMATIC                     1111
#     i = 0
#     species_lst = ['bulk','ooh','o','oh','bulk']
#     free_e_lst = []
#     for species in species_lst:
#         free_e_lst.append(fe_all[metal][species])
#     # print(free_e_lst)
#     #__|                                                                  1111
#
#     #| - Normalized energies by bulk energy                               1111
#     bulk_free_e =  free_e_lst[0]
#     for i in range(len(free_e_lst)):
#         free_e_lst[i] = free_e_lst[i] - bulk_free_e
#     #__|                                                                  1111
#
#     #| - Replaces String "M" with Metal Name                              1111
#     for j in range(len(species_lst)):
#         species_lst[j] = str(species_lst[j]).replace('M', metal)
#     #__|                                                                  1111
#
#     #| - Constructing Free Energy Diagram Data                            1111
#     xdata = []
#     ydata = []
#     for step in range(len(steps)):
#         xdata.append(steps[step])
#         xdata.append(steps[step]+1)
#
#         if steps[step] == 0:
#             ydata.append(free_e_lst[step] + 4.92)
#             ydata.append(free_e_lst[step] + 4.92)
#         else:
#             ydata.append(free_e_lst[step])
#             ydata.append(free_e_lst[step])
#     #__|                                                                  1111
#
#     #| - Finding Overpotential                                            1111
#     overpotential = 0
#     for i in range(len(ydata)-1):
#         if ydata[i+1] == ydata[i]:
#             None
#         else:
#             d = 1.23 + ydata[i+1] - ydata[i]
#             if d > overpotential:
#                 overpotential = d
#     #__|                                                                  1111
#
#     #| - Constructing FED at Equilibrium and Overpotential                1111
#     ydatalim    = []
#     ydataeq        = []
#     for j in range(len(free_e_lst)):
#         if j == 0:
#
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23 + 4*1.23)
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23 + 4*1.23)
#
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential)+4*1.23)
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential)+4*1.23)
#
#         else:
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23)
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23)
#
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential))
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential))
#
#     #__|                                                                  1111
#
#     #| - Constructing FED for H2O2 (2-e process)                          1111
#     ydata_h2o2    = []
#     for j in range(3):
#         if j == 0:
#             ydata_h2o2.append(free_e_lst[j] - (4-j)*(0.7/2.) + 4*1.23)
#             ydata_h2o2.append(free_e_lst[j] - (4-j)*(0.7/2.) + 4*1.23)
#
#         elif j == 2:
#             ydata_h2o2.append(3.52)
#             ydata_h2o2.append(3.52)
#
#         else:
#             ydata_h2o2.append(free_e_lst[j] - (4-2.*j)*(0.7/2.))
#             ydata_h2o2.append(free_e_lst[j] - (4-2.*j)*(0.7/2.))
#
#     xdata_h2o2 = [ 0 , 1 , 1 , 2 , 2 , 3 ]
#     #__|                                                                  1111
#
#     return {"ydata_nobias":ydata, "ydata_eq":ydataeq, "ydata_h2o2":ydata_h2o2,
#      "xdata":xdata,    "overpotential":overpotential, "species_fe":species_fe,
#     "xdata_h2o2":xdata_h2o2}
#
#     #__|
#
#
# def create_rxn_coord_array(rxn_steps, spacing=0, step_size=1):
#     """
#     Creates a reaction coordinate array ([0, 1, 1, 2, 2, 3]) for plotting
#
#     Args:
#         rxn_steps:
#             # <type 'int'>
#             # Number of steps in reaction coordinate including initial
#             and final state.
#             # Ex. A -> Intermediate -> C has 3 steps/states
#
#         spacing:
#             # <type 'float'>
#             # Spacing inbetween the energy levels. The default of 0 creates a
#             free energy diagram that looks like steps
#     """
#     #| - create_rxn_coord_array
#     lst = []
#     for i in range(1, rxn_steps):
#         if i == 1:
#             lst.append(step_size)
#             lst.append(step_size+spacing)
#         if i != 1:
#             lst.append(lst[-1]+step_size)
#             lst.append(lst[-2]+step_size+spacing)
#
#     lst.insert(0,0)
#     lst.append(lst[-1]+step_size)
#     return lst
#     #__|
#
#
# def plot_fed_line(plot_obj, *args, **kwargs):
#     """
#
#     Args:
#         plot_obj:
#             # <>
#             #
#         xdata:
#             # <type 'list'>
#             #
#         ydata:
#             # <type 'list'>
#             #
#     """
#     #| - plot_fed_line
#     xdata = args[0]; ydata = args[1]
#
#     plot_obj.plot(xdata[2*0:2*0+2], ydata[2*0:2*0+2], linewidth=3,**kwargs)
#
#     if "label" in kwargs: del kwargs["label"]
#
#     kwargs_even        = copy.copy(kwargs)
#     kwargs_odd        = copy.copy(kwargs)
#
#     if "path_effects" in kwargs_odd: del kwargs_odd["path_effects"]
#
#     for i in range(len(xdata)/2):
#         tmp = plot_obj.plot(xdata[2*i:2*i+2], ydata[2*i:2*i+2],
#         linewidth=3, solid_capstyle='round', **kwargs_even)
#         plot_obj.plot(xdata[2*i+1:2*i+3], ydata[2*i+1:2*i+3],
#         linewidth=0.5,**kwargs_odd)
#         # print(dir(tmp))
#         # tmp.set_solid_capstyle('round')
#         # ax.margins(.2)
#     #__|

#__|
