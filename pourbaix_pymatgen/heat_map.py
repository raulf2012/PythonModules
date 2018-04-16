#| -  - Import Modules
# -*- coding: utf-8 -*-
from pourdiag import pd_entries
from entry_methods import pure_atoms_remove, alloy_entries, base_atom
from pd_screen_tools import phase_coord, phase_filter
from stability_crit import most_stable_phase
from element_list import ElementList, ElemList_mod

import numpy as np
import matplotlib.colors as colors
#__|

def process(i,j,comp=0.5, heat_map_scale='Pt_ref', pt_oxid_V=0.6470339):
	"""
	Creates and analyzes a Pourbaix Diagram from two elements (they can be the
	same ex. Pt,Pt). Finds relevant entries, creates Pourbaix diagram,
	identifies stable phases, and calculates the stability of the phase

	Args:
		i: First element in the system
		j: Second element in the system
		comp: Composition loading of the two elements (default is 0.5)
	"""
	#| -  - process
	entries = pd_entries(i.symbol,j.symbol)
	coord = phase_coord(entries, comp, prim_elem=i.symbol)
	filt1 = phase_filter(coord,'metallic')
	filt2 = phase_filter(coord,'metallic_metallic')
	filt = filt1 + filt2

	# msp = most_stable_phase(filt,Pt_ref=True)
	if not filt:
		print 'heat_map.process - no phase present - '+i.symbol+'-'+j.symbol
		msp = [-1.5,'Pourbaix Entry placeholder']	# TEMP
	else:
		msp = most_stable_phase(filt, scale=heat_map_scale, pt_oxid_V=pt_oxid_V)
		# msp = most_stable_phase(filt,pH=10.5,scale=heat_map_scale)
	return msp
	#__|

def process_alloy(i,j):	# Deprecated *******************************************
	"""
	Creates and analyzes a Pourbaix Diagram from two elements (they can be the
	same ex. Pt,Pt). Finds relevant entries, removes the pure element entries,
	for each alloy removes all other alloys and analyzes stability of all
	forced alloy phases. Returns the the most highest performing forced alloy
	phase.

	Args:
		i: First element in the system
		j: Second element in the system
	"""
	#| -  - process_alloy
	entries = pd_entries(i.symbol,j.symbol)
#	entries = pure_atoms_remove(entries)
	alloy_entr = alloy_entries(entries)
	if not alloy_entr:
		print 'heat map - no alloy entries'
		non_alloy = process(i,j)
		return non_alloy

	base_atoms = base_atom(entries)
	alloy_performance_lst = []
	for alloy in alloy_entr:
		entries_0 = entries[:]
		comp = alloy.composition \
			.fractional_composition.get_atomic_fraction(base_atoms[0])
		for alloy_0 in alloy_entr:
			if not alloy_0 == alloy: entries_0.remove(alloy_0)

		coord = phase_coord(entries_0,comp)
		filt1 = phase_filter(coord,'metallic')
		filt2 = phase_filter(coord,'metallic_metallic')
		filt = filt1 + filt2
		try:
			alloy_performance_lst.append(most_stable_phase(filt, scale='Pt_ref'))
		except:
			pass

	alloy_performance_lst.sort(key=lambda x: x[0], reverse=True)
	# NOTE This sort will not work when the performance criteria is distance
	# from the ORR line (which needs to be as small as possible)
	try:
		best_alloy = alloy_performance_lst[0]
	except:
		best_alloy = [-1,'placeholder'] #NOTE Make this better
	return best_alloy
	#__|

def construct_output_matrix(elem):
	"""
	Constructs the skeleton of the output matrix from the element sweep. Matrix
	is square with dimensions of n, where n is the number of elements in elem

	Args:
		elem: List of elements to be screened over
	"""
	#| -  - construct_output_matrix

	o_lst = []
	i_cnt = 0
	for i in elem:
		o_lst.append([])
		for j in elem:
			o_lst[i_cnt].append([])
		i_cnt = i_cnt+1
	return o_lst
	#__|

def finish_symmetric_matrix(elem, output_lst):
	"""
	Fills in the the opposite diagonal of a symmetric matrix

	Args:
		elem: Element list being screened over
		output_lst: Half filled matrix to be filled in
	"""
	#| -  - finish_symmetric_matrix
	i_cnt = 0
	for i in elem[:-1]:
		j_cnt = 1
		for j in elem[1:]:
			output_lst[i_cnt][j_cnt] = output_lst[j_cnt][i_cnt]
			j_cnt = j_cnt+1
		i_cnt=i_cnt+1
	return output_lst
	#__|


################################################################################
#			██       ██████   ██████  ██████  ██ ███    ██  ██████
#			██      ██    ██ ██    ██ ██   ██ ██ ████   ██ ██
#			██      ██    ██ ██    ██ ██████  ██ ██ ██  ██ ██   ███
#			██      ██    ██ ██    ██ ██      ██ ██  ██ ██ ██    ██
#			███████  ██████   ██████  ██      ██ ██   ████  ██████
################################################################################

def run_all_binary_combinations(elements, loop_funct, scale, pt_oxid_V=0.6470339):
	"""
	"""
	#| -  - run_all_binary_combinations
	o_lst = construct_output_matrix(elements)

	print 'constructing output V_crit matrix from scratch'

	i_cnt = 0
	for j in elements:
	# for j in elements[1:]:	#TEMP
		# for i in elements[i_cnt:1]:	#TEMP
		for i in elements[i_cnt:]:
			print '##_'+j.symbol + i.symbol+'_##'
			output = loop_funct(i, j, heat_map_scale=scale, pt_oxid_V=pt_oxid_V)
			o_lst[elements.index(i)][elements.index(j)] = output

			#TEMP_PRINT
			if i.name=='Ru' or j.name=='Ru':
				print 'heat_map.run_all_binary_combinations - TEMP_PRINT'
				print output
			#

			print '_'
		print '_________________________________________'
		i_cnt=i_cnt+1
	print '#########################################'

	o_lst = finish_symmetric_matrix(elements,o_lst)

	return o_lst
	#__|

def oxidation_dissolution_product_0(i, j, scale):
	"""
	Creates Pourbaix Diagrams for single or binary systems
	"""
	#| -  - oxidation_dissolution_product_0
	# from pourdiag import pd_entries

	from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
	from pymatgen.analysis.pourbaix.plotter import PourbaixPlotter
	from pd_screen_tools import phase_coord, phase_filter
	from stability_crit import most_stable_phase, oxidation_dissolution_product

	elem0 = i.symbol; elem1 = j.symbol

	mat_co_0 = 0.50		# Composition of 1st entry in elem_sys

	entr = pd_entries(elem0,elem1)
	pourbaix = PourbaixDiagram(entr,{elem0: mat_co_0,elem1: 1-mat_co_0})
	plotter = PourbaixPlotter(pourbaix)

	coord = phase_coord(entr,mat_co_0)
	filt1 = phase_filter(coord,'metallic')
	filt2 = phase_filter(coord,'metallic_metallic')
	filt = filt1 + filt2
	msp = most_stable_phase(filt,scale='RHE')

	tmp = oxidation_dissolution_product(coord,msp)

	if 'Ion' in tmp:
		entry_lst = 'dis'
	else:
		entry_lst = 'oxi'

	"""
	entry_lst = ''
	i_cnt = 0
	for i in tmp:
		entry_lst = str(entry_lst)+'\n'+str(i)
		if i_cnt==0:
			entry_lst = entry_lst[1:]
		i_cnt=i_cnt+1

	"""
	return entry_lst
	#__|

# TEMP
def ref_atoms(i,j, scale):
	"""
	TEMP
	"""
	#| -  - ref_atoms
	ref_atoms = pd_entries(i.symbol,j.symbol)
	print ref_atoms	#TEMP_PRINT
	return ref_atoms
	print ref_atoms # TEMP_PRINT
	#__|

################################################################################
# 				██████  ██       ██████  ████████
# 				██   ██ ██      ██    ██    ██
# 				██████  ██      ██    ██    ██
# 				██      ██      ██    ██    ██
# 				██      ███████  ██████     ██
################################################################################

class MidpointNormalize(colors.Normalize):
	"""
	Used with diverging color schemes to set the white color to 0
	"""
	#| -  - MidpointNormalize
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))
	#__|

def extract_data_from_matrix(output_lst):
	"""
	Extracts numerical data from screen output matrix (because the output
	matrix entries contain numerical stability and PD entry data) and makes a
	completley numerical matrix for plotting

	Args:
		output_lst: Matrix containing numerical data paired with entry data
	"""
	#| -  - extract_data_from_matrix
	data_matrix=[]
	cnt=0
	for i in output_lst:
		data_matrix.append([])
		for j in i:
			try:
				data_matrix[cnt].append(j[0])
			except:
				data_matrix[cnt].append(0)
		cnt=cnt+1
	return data_matrix
	#__|

def plot_heat_map(data_matrix, elem, text_overlay=None, composition=False,
	show_plot=False, save_file=True, file_type='.pdf',
	heat_map_scale='Pt_ref', transparency=False, plot_title=None,
	colorbar_title=None, lock_cbar_rnge=None):
	"""
	Constructs heat map plot from a matrix. If another matrix is passed to
	text_overlay it will be overlayed as text on top of the heatmap squares

	Args:
		data_matrix:
		elem:
		text_overlay: Matrix of data which will be overlayed on the heatmap, if
			='data_value' it will overlay the numerical value for each grid point
	"""
	#| -  - plot_heat_map

	#| -  - Import modules
	from heat_map import MidpointNormalize
	import matplotlib.pyplot as plt; import numpy as np
	from element_list import elem_str_mke
	#__|

	#| -  - Setting colorbar range
	if lock_cbar_rnge == None:
		vmin_=None
		vmax_=None
	elif lock_cbar_rnge != None:
		vmin_=lock_cbar_rnge[0]
		vmax_=lock_cbar_rnge[1]
	#__|

	#| -  - Main, create heatmap, set labels, colorbar, scale plot
	font_color = 'White'

	fig, ax1 = plt.subplots(1,1)

	if heat_map_scale=='Pt_ref':
		cmap_ = 'seismic'
		norm_ = MidpointNormalize(midpoint=0.)
	elif heat_map_scale=='RHE':
		cmap_ = 'hot'
		norm_ = None

	img = ax1.imshow(	data_matrix, cmap=cmap_, interpolation='nearest',
						norm=norm_, vmin=vmin_, vmax=vmax_)

	#| -  - TEMP - Plotting heat map with either Pt_ref or RHE scaling
	# if heat_map_scale=='Pt_ref':
	# 	cmap_='seismic'
	# 	img = ax1.imshow(	data_matrix, cmap=cmap_, interpolation='nearest',
	# 					norm=MidpointNormalize(midpoint=0.), vmin=vmin_, vmax=vmax_)
	#
	# elif heat_map_scale=='RHE':
	# 	cmap_='seismic'
	# 	# cmap_='hot'
	#
	# 	img = ax1.imshow(	data_matrix, cmap=cmap_, interpolation='nearest',
	# 					norm=MidpointNormalize(midpoint=0.), vmin=vmin_, vmax=vmax_)

		# img = ax1.imshow(	data_matrix, cmap=cmap_, interpolation='nearest',
		#  				norm=MidpointNormalize(midpoint=0.))					#TEMP
		# img = ax1.imshow(data_matrix, cmap=cmap_, interpolation='nearest')
	#cmap='hot'; cmap='seismic'; cmap='Reds'
	#__|


	elem_str = elem_str_mke(elem)
	ax1.tick_params(labelbottom='off',labeltop='on')
	# plt.setp(ax1.get_xticklabels(), fontdict=font)
	plt.setp(ax1.get_xticklabels(), fontsize=32, color=font_color)
	plt.setp(ax1.get_yticklabels(), fontsize=32, color=font_color)
	tcks = range(0,len(elem_str))
	plt.xticks(tcks,elem_str)
	plt.yticks(tcks,elem_str)

	if plot_title != None:
		plt.title(plot_title, y=1.06, fontsize=25./23*len(elem)+13,
					color=font_color)

	# fig_note = 	'Note: Text box numbers indicate the number of alloys
	# available in the Materials\nProject data base for a given binary
	# element pair. Black boxes indicate no alloy.'

	# plt.figtext(0.08,0.08, fig_note,fontsize=26)

############################ COLORBAR ##########################################
	cb = plt.colorbar(img, spacing='uniform', fraction=0.046, pad=0.04)
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=font_color)
	# cb.ax.tick_params(labelsize=20./23*len(elem)+5)
	cb.ax.tick_params(labelsize=20./23*len(elem)+14)

	## Setting the colorbar label
	if colorbar_title != None:
		col_bar_lb = colorbar_title
	elif heat_map_scale=='Pt_ref':
		col_bar_lb = 'V$\mathregular{_{critical}}$ - V$\mathregular{_{Pt}}$ [V vs RHE]'
	elif heat_map_scale=='RHE':
		col_bar_lb = 'V$\mathregular{_{critical}}$ [V vs RHE]'

	cb.set_label(col_bar_lb, rotation=-90, fontsize=25./23*len(elem)+10,
	labelpad=50, color=font_color)
################################################################################

	scl = 20./23.*len(elem)+2; fig.set_figheight(scl); fig.set_figwidth(scl)
	#__|

	if composition != False:
		# Adding text to plot
		import matplotlib.patheffects as PathEffects

		plt.text(0.1, 0.9,str(composition), fontsize=45, ha='center',
		va='center', color=font_color,
		transform=ax1.transAxes).set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

	#| -  - Setting text in heatmaps squares
	###################### SETTING TEXT IN EACH SQUARE #########################
	############################################################################
	if text_overlay == 'data_value':
		text_overlay = np.array(data_matrix)

		diff = 1.
		min_val = 0.
		rows = text_overlay.shape[0]
		cols = text_overlay.shape[1]

		col_array = np.arange(min_val, cols, diff)
		row_array = np.arange(min_val, rows, diff)
		x, y = np.meshgrid(col_array, row_array)
		# TEMP
		print '_________________________________________'
		#
		import matplotlib.patheffects as PathEffects
		for col_val, row_val in zip(x.flatten(), y.flatten()):
			c = np.round(text_overlay[int(row_val), int(col_val)],2)
			fontsize_0=17

			ax1.text(col_val, row_val, c, fontsize=fontsize_0 ,va='center',
			ha='center').set_path_effects([PathEffects.withStroke(linewidth=2,
			foreground='w')])

	elif text_overlay != None:
		text_overlay = np.array(text_overlay)

		diff = 1.
		min_val = 0.
		rows = text_overlay.shape[0]
		cols = text_overlay.shape[1]

		col_array = np.arange(min_val, cols, diff)
		row_array = np.arange(min_val, rows, diff)
		x, y = np.meshgrid(col_array, row_array)

		import matplotlib.patheffects as PathEffects
		for col_val, row_val in zip(x.flatten(), y.flatten()):

			c = str(text_overlay[row_val.astype(int),col_val.astype(int)])
			fontsize_0=54/2

			#| -  - TEMP
			# if row_val.astype(int)==col_val.astype(int):
			# 	c = ''
			# 	fontsize_0=10
			# elif text_overlay[row_val.astype(int),col_val.astype(int)] == 0:
			# 	c = u'\u25A0'
			# 	fontsize_0=20 # 56 covers most of square
			# else:
			# 	c = str(text_overlay[row_val.astype(int),col_val.astype(int)])
			# 	fontsize_0=26/4
			#__|

			ax1.text(col_val, row_val, c, fontsize=fontsize_0,
			color=font_color, va='center',
			ha='center').set_path_effects([PathEffects.withStroke(linewidth=2,
			foreground='black')])
	############################################################################
	#__|

	#| -  - Saving Figure, showing figure
	import fnmatch; import sys; import os
	num_fle = len(fnmatch.filter(os.listdir('.'),'*'+file_type))
	fle_nme = 'fig_heat_map'+'_'+str(num_fle)+file_type

	if save_file == True:
		print 'saving '+fle_nme

		# if file_type=='.svg':
		# 	fig.savefig(fle_nme, format='svg',spi=1200,transparent=transparency)
		# else:
		# 	fig.savefig(fle_nme, transparent=transparency)

		fig.savefig(fle_nme, format=file_type[1:],spi=1200,transparent=transparency)

	if show_plot==True: fig.patch.set_facecolor('black'); plt.show()
	#__|

	#__|
