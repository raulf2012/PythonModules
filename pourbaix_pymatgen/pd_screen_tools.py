def ORR_line(pH):
	"""
	"""
	#| -  - ORR_line
	intercept = 1.23
	slope = -0.0591

	V = slope*pH + intercept
	return V
	#__|

def stable_mat_one_elem(pH,V,mat):
	"""
	Returns pymatgen pourbaix entry object corresponding to most stable species

	Args:
		pH: pH of system
		V: Potential with reference to the Standard Hydrogen Electrode (SHE)
	"""
	#| -  - stable_mat_one_elem
	from pourdiag_single import pourdiag_single as pd_sgle
	from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
	# Access entry Gibbs(pH, V)
	from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer


	pH = 0
	V = 1
	mat = 'Pt'

	all_entries = pd_sgle(mat)

	pourbaix = PourbaixDiagram(all_entries)
	PA=PourbaixAnalyzer(pourbaix)

	templist=[]
	for i in all_entries:
		templist.append(PA.g(i,pH,V))
	minE=min(templist)

	for i in all_entries:
		if PA.g(i,pH,V)==minE:
			StableSpec=i
	return StableSpec		# Returns the entries object of the stable species

	#__|

def plot_reg(coord_data):
	"""
	Plots the region boundaries for the given region

	Args:
		coord_data: Coordinate data of region of interest, must be in the form
		of [[[x-points],[y-points]],[], ...]
	"""
	#| -  - plot_reg
	import numpy as np
	import matplotlib.pyplot as plt

	fig, ax = plt.subplots()
	for line in coord_data:
		ax.plot(line[0],line[1])
	plt.show()
	#__|

def phase_coord(entries, atom_comp, prim_elem=False):
	"""
	Produces a list of line segments corresponding to each phase area in a PD
	along with the PD entries corresponding to each area.
	The produced list is of the following form:
		list = [[[coordinate data], [pourbaix entry]], [], [] .... ]

	Args:
		entries: List of entries in a PD
		atom_comp: Composition of atom if system is binary, given as a fraction
			between 0 and 1. Corresponds to the element with lowest atomic number if
			prim_elem is left to its default
		prim_elem: Primary element to which the atom_comp is assigned
	"""
	#| -  - phase_coord
	from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
	from pymatgen.analysis.pourbaix.plotter import PourbaixPlotter
	from entry_methods import base_atom

	base_atoms = base_atom(entries)
	mat0 = base_atoms[0]
	if len(base_atoms)==2: mat1 = base_atoms[1]
	else: mat1 = mat0

	# If the primary element is declared then set the given composition to it
	if not prim_elem==False:
		for atom in base_atoms:
			if atom==prim_elem:
				mat0 = atom
			else: mat1 = atom

	pd = PourbaixDiagram(entries,{mat0: atom_comp,mat1: 1-atom_comp})
	pl = PourbaixPlotter(pd)
	ppd = pl.pourbaix_plot_data([[-2, 16],[-3, 3]])	#ppd == Pourbaix_Plot Data
	pd_lst = []
	cnt = 0

	for stable_entry in ppd[0]:
		pd_lst.append([])
		pd_lst[cnt].append(ppd[0][stable_entry])
		pd_lst[cnt].append(stable_entry.entrylist)
		cnt = cnt+1

	return pd_lst

	#__|

def phase_filter(phase_coord,phase_type):
	"""
	Returns a list of Pourbaix diagrams regions and corresponding species
	that match the specified phase_type

	Args:
		phase_coord: PD phase coordinate data produced from phase_coord
		phase_type: Type of phase that will be filtered for. Samples include the following:
				metallic, oxide, metallic_metallic, metallic_oxide, oxide_oxide,
				metallic_aqueous oxide_aqueous, aqueous_aqueous
	"""

	#| -  - phase_filter
	## METALLIC: 1 METALLIC PHASE
	#	For a binary system this is necessarily an alloy
	#TMP

	if phase_type == 'metallic':
		met_phase_lst = []
		for region in phase_coord:
			is_oxide_phase = False
			if len(region[1]) == 1:
				if region[1][0].phase_type == 'Solid':
					for elem in region[1][0].composition.elements:
						if elem.symbol == 'O': is_oxide_phase = True
					if is_oxide_phase == False:
						met_phase_lst.append(region)
		return met_phase_lst

	## METALLIC_METALLIC: 2 METALLIC SPECIES
	#	May be 2 single atomic species (ex. Ni(s) + Pt(s)) or two alloys (ex. NiPt2(s) + Ni2Pt(s))

	if phase_type == 'metallic_metallic':
		met_met_phase_lst = []

		for region in phase_coord:
			c1 = len(region[1]) == 2

			if len(region[1]) == 2:
				c2 = region[1][0].phase_type == 'Solid'
				c3 = region[1][1].phase_type == 'Solid'

				if c2 and c3:
					is_oxide_phase = False

					for elem in region[1][0].composition.elements:
						if elem.symbol == 'O': is_oxide_phase = True
					for elem in region[1][1].composition.elements:
						if elem.symbol == 'O': is_oxide_phase = True

					if is_oxide_phase == False:
						met_met_phase_lst.append(region)
		return met_met_phase_lst

	#__|

def is_solid_phase(mat1, mat2, mat1_co=0.5):
	"""
	Returns TRUE is there exists a all solid phase in the binary Pourbaix Diagram
	This means that the phase doesn't have any aqueous species

	Args:
		mat1:
		mat2:
		mat1_co:
	"""

	#| -  - is_solid_phase
	from pourdiag import pourdiag # Returns Pourbaix entries for binary system
	from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
	from pymatgen.analysis.pourbaix.plotter import PourbaixPlotter

	mat2_co = 1-mat1_co

	pd_b = pourdiag(mat1,mat2)
	return pd_b

	pd = PourbaixDiagram(pd_b,{mat1: mat1_co,mat2: mat2_co})
	pl = PourbaixPlotter(pd)
	ppd = pl.pourbaix_plot_data([[-2, 16],[-3, 3]])		#ppd == Pourbaix_Plot Data

	pd_lst = []
	cnt = 0
	for stable_entry in ppd[0]:
		pd_lst.append([])
		pd_lst[cnt].append(ppd[0][stable_entry])
		pd_lst[cnt].append(stable_entry.entrylist)
		cnt = cnt+1

	solidphase=False
	for i in pd_lst:
		if len(i[1])==1:
			if i[1][0].phase_type == 'Solid':
				solidphase=True
		if len(i[1])==2:
			if i[1][0].phase_type and i[1][1].phase_type == 'Solid':
				solidphase=True
	return solidphase

	#__|
