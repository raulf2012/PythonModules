def ref_atoms_dict():
	"""
	Contains dictionary of reference atom energies and entropic corrections to
	form Gibbs free energies

	Args:

	"""
	#start_fold - ref_atoms_dict

	# e_o = -4.93552791875		# mp-12957, half of O2 - No entropy term
	# e_h = -3.2397				# mp-754417, half of H2 - No entropy term

	# e_o = -5.25225891875		# mp-12957, half of O2 - With entropy @300K
	# e_h = -3.6018845			# mp-754417, half of H2 - With entropy @300K

	# e_o = -5.11					# Optimzed for 1st row transitino metals
	# e_h = -3.6018845			# mp-754417, half of H2 - With entropy @300K

	e_o = -5.21
	e_h = -3.6018845

	#start_fold - H2O, H2 Reference State
	# From Michal Badjich
	o2=-9.88557216
	h2o=-14.23091949
	h2=-6.77149190
	h2ozpe=0.5584
	h2ocp=0.10
	h2zpe=0.27283
	h2cp=0.09
	o2zpe=0.098
	o2cp=0.090

	# o_ref1=o2/2
	o_ref2=(h2o-h2 +2.506)+(h2ozpe-h2zpe-0.5*o2zpe)+(h2ocp-h2cp-0.5*o2cp)

	# e_o = o_ref2
	# e_h = h2/2.
	#end_fold


	#start_fold - Entropic Corrections
	entr_o = -0.3167
	entr_h = -0.36218


	ref_dict = {}
	ref_dict['e_o'] = e_o
	ref_dict['e_h'] = e_h

	return ref_dict

	#end_fold

	#end_fold

def h2o_energy():
	"""
	Enthalpy and Gibbs free energy of formation for liquid water
	"""
	e_h2o = -2.96211
	g_h2o = -2.4583

	h2o_dict = {}
	h2o_dict['h2o_enth'] = e_h2o
	h2o_dict['h2o_gibb'] = g_h2o

	return h2o_dict

def form_e_all_oxides(element, only_pd_entries=False):
	# e_o = oxy_energy
	"""
		OUTLINE:
		1. Get all entries in the chemical ensemble for 1 element and oxygen
		2. Extract all of the oxide species
		3. Extract the pure metallic references
		4. Apply correction scheme to entries
		5. Calculate the formation energy for the oxides per oxygen molecule
	"""
	#start_fold - form_e_all_oxides

	#start_fold - Imported Modules
	from pd_make import entry_data, aq_correction, stable_entr, form_e
	from entry_methods import entry_remove, pure_atoms_return,contains_element, norm_e
	from energy_scheme import ref_atoms_dict
	#end_fold

	#start_fold - Script parameters
	mprester_key = 'ZJhfHmMTTwbW29Sr'	# Input your materials project id
	direct_0 = '/home/flores12/01_ORR-MatStabScreen/01_virenv-pymatgen/01_data/01-1_local_MP_entry/'
	correction_=True
	#end_fold

	all_entries = entry_data(element,element,direct_0,mprester_key)['entries']
	oxygen_entries = contains_element(all_entries,'O')
	hydrogen_entries = contains_element(all_entries,'H')

	entries_no_h = [entry for entry in oxygen_entries if entry not in hydrogen_entries ]

	oxides = entry_remove(entries_no_h, 'O2')

	#start_fold - Finding the Entries Used in the Pourbaix Diagram
	entry_ion_data = entry_data(element, element, direct_0, mprester_key)
	entries = entry_ion_data["entries"]

	entries_aqcorr = aq_correction(entries)
	stable_solids_minus_h2o = stable_entr(entries_aqcorr)
	pbx_solid_entries = form_e(stable_solids_minus_h2o, entries_aqcorr)

	entry_id_lst = []
	for i in pbx_solid_entries:
		entry_id_lst.append(i.entry_id)
	#end_fold


	#start_fold - Deleting Entries in Oxides Which Aren't in Pourbaix Entry List
	if only_pd_entries==True:
		oxides_temp = []
		for i in oxides:
			if i.entry_id in entry_id_lst:
				oxides_temp.append(i)
		oxides = oxides_temp
	#end_fold

	#start_fold - Metallic References
	elem_ref = pure_atoms_return(all_entries)[0]
	elem_ref_e = norm_e(elem_ref)

	if not elem_ref.name==element:
		print 'calc_form_e - The reference atom is not the same as the element'
	#end_fold

	#start_fold - Formation Energy
	ref_atom_energies = ref_atoms_dict()
	e_o = ref_atom_energies['e_o']
	e_h = ref_atom_energies['e_h']

	form_e_lst = []
	for oxide in oxides:
		oxide_e = norm_e(oxide,correction=correction_)

		elem_coef = oxide.composition.get_el_amt_dict()[element]/oxide.composition.get_integer_formula_and_factor()[1]
		oxygen_coef = oxide.composition.get_el_amt_dict()['O']/oxide.composition.get_integer_formula_and_factor()[1]

		form_e = (oxide_e - elem_coef*elem_ref_e - oxygen_coef*e_o)/(oxygen_coef/2.)

		elem_comp_lst = []
		elem_comp_lst.append(str(element)+str(elem_coef)[:1])
		elem_comp_lst.append('O'+str(oxygen_coef)[:1])

		entry_dict = {}
		entry_dict['form_e'] = form_e
		entry_dict['entry'] = oxide.name
		entry_dict['entry_elements'] = elem_comp_lst

		form_e_lst.append(entry_dict)

	#end_fold
	return form_e_lst

	#end_fold
