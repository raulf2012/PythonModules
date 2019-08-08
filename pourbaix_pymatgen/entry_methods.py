# -*- coding: utf-8 -*-

# ███    ███ ██████      ██████   █████  ████████ ██████  ███████ ███████
# ████  ████ ██   ██     ██   ██ ██   ██    ██    ██   ██ ██      ██
# ██ ████ ██ ██████      ██   ██ ███████    ██    ██████  ███████ █████
# ██  ██  ██ ██          ██   ██ ██   ██    ██    ██   ██      ██ ██
# ██      ██ ██          ██████  ██   ██    ██    ██████  ███████ ███████


def get_entries_MP(material_name):
	"""
	Args:
		material_name: Material's formula unit
	"""
	#| -  - get_entries_MP
	from pymatgen.matproj.rest import MPRester # MP API key to access species in MP
	import warnings
	warnings.filterwarnings('ignore')
	mpr = MPRester('ZJhfHmMTTwbW29Sr')

	# entries = mpr.get_entries
	#__|

def get_entry_MP(entry_s):
	"""
	Returns a specific entry from the Materials Project database
	Args:
		ent:
	"""
	#| -  - get_entry_MP
	from pymatgen.matproj.rest import MPRester # MP API key to access species in MP
	import warnings
	warnings.filterwarnings('ignore')
	mpr = MPRester('ZJhfHmMTTwbW29Sr')

	# entry = mpr.get_entries(entry_s)[0]


	entryid_prefix = entry_s[:3]
	if entryid_prefix=='mp-':
		out = mpr.get_entries(entry_s)[0]

	else:
		out = mpr.get_entries(entry_s)
	return out
	#__|

def norm_e(entry, correction=True):
	"""

	Args:
		entry:
		correction: Apply default MP corrections
	"""
	#| -  - norm_e
	if correction==True:
		e_raw = entry.energy
	elif correction==False:
		e_raw = entry.uncorrected_energy

	stoich_fact = entry.composition.get_integer_formula_and_factor()[1]

	e_norm = e_raw/stoich_fact

	return e_norm
	#__|

def return_entry(entry_list, entryid):
	"""
	Returns the desired entry from a list of pymatgen entries if entryid given
	Returns a list of entries matching the given entry name
	Args:
		entry_list: List of pymatgen entries
		entryid: Materials Project id (ex. mp-715572) or entry name
	"""
	#| -  - return_entry
	entryid_prefix = entryid[:3]
	if entryid_prefix=='mp-':
		for ent in entry_list:
			if ent.entry_id == entryid:
				return ent
	else:
		out_lst = []
		for ent in entry_list:
			if ent.name == entryid:
				out_lst.append(ent)
		return out_lst

	#__|


def entry_remove(entries, entry_to_remove):
	"""
	Removes the entries in a list which match the name of entry_to_remove

	Args:
		entries: List of entries to be processed
		entry_to_remove: Name of entry to be removed from the entries list
	"""
	#| -  - entry_remove
	import numpy as np

	index_lst = []
	entries_copy = entries[:]
	cnt = 0
	for entry in entries_copy:
		if entry.name == entry_to_remove:
			index_lst.append(cnt)
			# entries.remove(entry)
		cnt=cnt+1

	entries_copy = np.delete(entries_copy, index_lst).tolist()

	return entries_copy
	#__|

def base_atom(entries):
	"""
	Returns the base metal(s) from a list of entries as a list of strings in
	order of their atomic number

	Args:
	    entries: List of entries
	"""
	#| -  - base_atom
	from pymatgen import Element

	elem_lst = []
	# Adds elements from every element to a list (losts of duplicates)
	for entry in entries:
		elem_lst.extend(entry.composition.elements)
	elem_lst_duplrem = list(set(elem_lst))	# Remove duplicate elements
	elem_lst_duplrem.sort(key=lambda x: x.number) # Order by atomic number

	# Converts list of element objects to atomic symbol strings
	elem_lst_nme = []
	for elem in elem_lst_duplrem:
		elem_lst_nme.append(elem.name)

	try:							# Attempts to remove oxygen from list
		elem_lst_nme.remove('O')	# if no oxygen is present do nothing
	except:
		pass
	try:							# Attempts to remove hydrogen from list
		elem_lst_nme.remove('H')	# if no hydrogen is present do nothing
	except:
		pass

	elem_lst = elem_lst_nme

	return elem_lst
	#__|

def pure_atoms_remove(entries):
	"""
	Removes entries corresponding to pure atoms, defined as entries with one
	species, and that species has to be the base element (or one of the base
	elements if binary)

	Args:
		entries: List of entries
	"""
	#| -  - pure_atoms_remove
	from pymatgen import Element
	from entry_methods import base_atom

	base_atoms = base_atom(entries)
	pure_met_lst=[]
	# Selects single atom entries of either elem_0 or elem_1 (not ion species)
	for entry in entries:
		entr_elem_lst = entry.composition.elements

		if 	len(entr_elem_lst)==1 and entry.phase_type!='Ion' \
				and entr_elem_lst[0].name in base_atoms:
			pure_met_lst.append(entry)
	# Removes the entries from entry list
	for entry in pure_met_lst:
		entries.remove(entry)
	return entries
	#__|

def pure_atoms_return(entries):
	"""
	"""
	#| -  - pure_atoms_return
	from pymatgen import Element
	from entry_methods import base_atom

	base_atoms = base_atom(entries)
	pure_met_lst=[]

	for entry in entries:
		entr_elem_lst = entry.composition.elements


		if 	len(entr_elem_lst)==1 and entr_elem_lst[0].name in base_atoms:
			pure_met_lst.append(entry)

	#| -  - If there are more pure atomic species than there are base atoms return an error
	if len(pure_met_lst)>len(base_atoms):
		tmp = 4
		# print 'entry_methods.pure_atoms_return - There is more than one option for the reference metallic atom'
	#__|

	# for i in pure_met_lst:
	# 	print norm_e(i)
	# 	print ''

	return pure_met_lst
	#__|


	# ████████ ███████ ███    ███ ██████
	#    ██    ██      ████  ████ ██   ██
	#    ██    █████   ██ ████ ██ ██████
	#    ██    ██      ██  ██  ██ ██
	#    ██    ███████ ██      ██ ██


def alloy_entries(entries):
	"""
	Returns the alloy species in entries in a list
	Args:
	    entries: List of entries
	"""
	#| -  - alloy_entries
	from pymatgen import Element
	from entry_methods import base_atom

	base_atoms = base_atom(entries)
	if len(base_atoms)==1:
		alloy_lst=[]
		return alloy_lst
	elem_0 = base_atoms[0]
	elem_1 = base_atoms[1]

	alloy_lst = []
	for entry in entries:
		entr_elem_lst = entry.composition.elements	# Elements in current entry
		entr_elem_lst_str = []
		for i in entr_elem_lst:
			entr_elem_lst_str.append(i.name)
		if len(entr_elem_lst_str)==2 and \
		(elem_0 in entr_elem_lst_str and elem_1 in entr_elem_lst_str):
			alloy_lst.append(entry)

	# Orders alloys from lowest to highest composition of elem_0
	alloy_lst.sort(key=lambda x: \
	x.composition.fractional_composition.get_atomic_fraction(elem_0))

	return alloy_lst
	#__|


def contains_element(entries,element_symbol):
	"""
	Returns the entries in a list of entries which contain element_symbol

	Args:
		entries: List of pymatgen objects
	"""
	#| -  - contains_element

	from pymatgen import Element


	#| -  - Get rid of entries without oxygen
	element = Element(element_symbol)

	lst = [entry for entry in entries if element in entry.composition.elements]
	# for i in entries:
	return lst
	#__|

	#__|


def number_of_elem_in_entry(entry,element):
	"""
	Returns the number of a certain element in an entry, where the entry has been
	normalized to a reduced formula unit

	Args:
		entry: entry/species to be processes
		element: element to be counted in the target entry
	"""
	#| -  - number_of_elem_in_entry
	#__|


def element_formula_list(entry):
	"""
	Returns the reduced formula list for an entry

	Args:
		entry: entry to be proccessed
	"""
	#| -  - element_formula_list
	formula_list = []
	for element in entry.composition.get_el_amt_dict():
		elem_coef = entry.composition.get_el_amt_dict()[element]/entry.composition.get_integer_formula_and_factor()[1]
		formula_list.append(str(element)+str(elem_coef)[:1])
	return formula_list
	#__|
