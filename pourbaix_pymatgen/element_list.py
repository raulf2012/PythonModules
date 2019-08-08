from pymatgen.core.periodic_table import Element	# Call elements by atomic #
class ElementList(object):
	"""
	Class for creating element lists and ranges
	By default the first list created is in order of atomic number
	Args:
	    N/A
	"""
	#| -  - ElementList
	def __init__(self, atom_num_lst=[[]]):
		self.atom_num_lst = atom_num_lst
		if not atom_num_lst==[[]]:
			self.list = self.mk_lst_atnum()
		self.trans_met = self.mk_lst_trans_met()

	def mk_lst_atnum(self):
		"""
		Makes list of elements from the atom_num_lst input
		"""
		elem_rnge=[]
		for i in self.atom_num_lst:
			el_strt=i[0]
			el_end=i[1]
			rnge_sect=range(el_strt,el_end+1)
			elem_rnge.extend(rnge_sect)
		elements=[]
		for i in elem_rnge:
			element=Element.from_Z(i)	# Indice -> pymatgen element object
			elements.append(element)
		return elements
		print elements

	def mk_lst_trans_met(self):
		"""
		Produces list of transition metals in order of atomic number
		"""
		elem_rnge_I = [[21,30],[39,44],[46,48],[74,76],[78,80]]
		elem_rnge=[]
		for i in elem_rnge_I:
			el_strt=i[0]
			el_end=i[1]
			rnge_sect=range(el_strt,el_end+1)
			elem_rnge.extend(rnge_sect)
		elements=[]
		for i in elem_rnge:
			element=Element.from_Z(i)	# Indice -> pymatgen element object
			elements.append(element)
		return elements

	#__|

class ElemList_mod(object):
	"""
	"""
	#| -  - ElemList_mod
	def __init__(self,elem_lst):
		self.elem_lst = elem_lst
	@property
	def sort(self,srt_type='group_num'):
		"""
		Sorts element entries
		Args:
		    srt_type
			DEFAULT: group_num, sorts by group number with lowest group number first, and
			lowest period first
		"""
		if srt_type == 'group_num':
			elem_lst = self.elem_lst
			elem_lst_sorted = elem_lst[:]		# Maintains the original unaffected
			elem_lst_sorted.sort(key=lambda x: x.group)
			return elem_lst_sorted

		#| -  - d-band filling
		if srt_type == 'd-band': #COMBAK
			tmp = 7




		#__|

	def remove(self,element):
		"""
		Removes element from element list
		Args:
			List of element objects
			If only 1 element just 1 element number or chemical symbol
			name "or" atomic number of the element to be removed
		"""
		elem_lst_0 = self.elem_lst
		if type(element)==type([]):
			for elem in element:
				if type(elem)==type('string'):
					# Find the occurance of input 'elem' in the list and returns it
					elem_0 = next((x for x in elem_lst_0 if x.name == elem), None)
					elem_lst_0 = [x for x in elem_lst_0 if not x.name == elem]
				elif type(elem)==type(2):
					elem_0 = next((x for x in elem_lst_0 if x.number == elem), None)
					elem_lst_0 = [x for x in elem_lst_0 if not x.number == elem]
			return elem_lst_0
		if type(element)==type('string'):
			# Find the occurance of input 'elem' in the list and returns it
			elem_0 = next((x for x in self.elem_lst if x.name == element), None)
			elem_lst_new = [x for x in self.elem_lst if not x.name == element]
			return elem_lst_new

		elif type(element)==type(2):
			elem_0 = next((x for x in self.elem_lst if x.number == element), None)
			elem_lst_new = [x for x in self.elem_lst if not x.number == element]
			return elem_lst_new
	#__|

def elem_str_mke(elem_lst):
	"""

	Args:
		elem_lst:
	"""
	#| -  - elem_str_mke
	elem_str=[]
	for i in elem_lst:
		elem_str.append(i.symbol)
	return elem_str

	#__|
