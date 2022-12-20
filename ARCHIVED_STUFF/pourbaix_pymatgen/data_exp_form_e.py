"""TEMP TEMP - 180317"""

def formation_e_data_dict():
	"""
	Dictionary of formation energies for transition metal species.

	All energies are in [eV]
	"""
	# | -  - formation_e_data_dict
	import ast
	direct_0 = '/home/flores12/01_ORR-MatStabScreen/01_virenv-pymatgen/01_data/03_formation_e/'
	List = open(direct_0+'data.py').read().splitlines()
	# Returns python list of entries_lst from read data file
	entries_lst = []
	line_cnt = 1
	for i in List:
		# print str(line_cnt)+'	:'+str(i)
		if not i=="":
			if i[0]=='#':
				line_cnt=line_cnt+1
				continue
			else:
				try:
					entries_lst.append(ast.literal_eval(i))
				except:
					print 'data_exp_form_e.formation_e_data_dict - error with line: '+str(line_cnt)
		line_cnt=line_cnt+1
	## Turns entry lists into dictionaries
	cf_k = "chem_formula"	# Chemical formula list key
	fe_k = "form_e"			# Formation energy key
	r_k = "reference"		# Literature reference key
	c_k = "comments"		# Comments key

	entries = []
	for i in entries_lst:
		entry_dict = {}
		entry_dict[cf_k] = i[0]
		entry_dict[fe_k] = i[1]
		entry_dict[r_k] = i[2]
		entry_dict[c_k] = i[3]
		entries.append(entry_dict)

	chem_form_lst = []
	for entry in entries:
		# chem_form = entry[0]
		chem_form_lst.append(frozenset(entry[cf_k]))
	unique_entries_lst = list(list(i) for i in list(set(chem_form_lst)))

	dict = {}
	for i in unique_entries_lst:
		dict[frozenset(i)]=[]
		for j in entries:
			if set(j[cf_k])==set(i):
				dict[frozenset(i)].append(j)
	return dict

	# __|

def get_entry(chemical_formula_list, data_dict):
	"""

	Args:
		chemical_formula_list:
		data_dict:
	"""
	# | -  - get_entry
	# print chemical_formula_list
	key = frozenset(chemical_formula_list)

	try:
		entry = data_dict[key]
	except KeyError:
		entry = None
		# print 'data_exp_form_e.get_entry - no experimental data available'

	return entry

	# __|

# form_e_data = formation_e_data_dict()
# temp = get_entry(['Pt1','O2'],form_e_data)

def oxygen_stoich(exp_entry):
	"""

	Args:
		exp_entry
	"""
	# | -  - oxygen_stoich
	formula = exp_entry['chem_formula']
	for i in formula:
		if i[0]=='O' and i[1].isdigit()==True:
			num_of_digits = len(i)-1
			if num_of_digits==1:
				oxy_num = int(i[1])
			elif num_of_digits>1:
				oxy_num = int(i[1:])
			break
	return oxy_num

	# __|


# | -  - Possible code to take "Pt2O3 and extract out the elements and stoicheometric numbers"
# import re
#
# chem_form_lst = []
# for entry in entries_lst:
#
# 	chem_form = entry[0]
# 	form_e = entry[1]
# 	ref = entry[2]
# 	comments = entry[3]
#
# 	chem_form_lst.append(set(chem_form))
#
# 	## Separating the element name from the stoicheometric factor
# 	for element in chem_form:
# 		elem_stoich = re.match(r"([a-z]+)([0-9]+)", element, re.I).groups()
# 		elem = elem_stoich[0]
# 		stoich = elem_stoich[1]
# __|
