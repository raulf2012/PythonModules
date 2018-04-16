# -*- coding: utf-8 -*-
from pd_make import entry_data, aq_correction, stable_entr, form_e, mke_pour_ion_entr

def pd_entries(mtnme_1,mtnme_2):
	"""
	Creates the entry objects corresponding to a binary or single component
	Pourbaix diagram

	Args:
		mtnme_1: Name of element 1
		mtnme_2: Name of element 2
	"""

	################################## INPUTS #######################
	mprester_key = 'ZJhfHmMTTwbW29Sr'	# Input your materials project id
	# Local directory containing entry data (solid and ion)
	direct_0 = '/home/flores12/01_ORR-MatStabScreen/01_virenv-pymatgen/01_data/01-1_local_MP_entry/'
	#################################################################
	entry_ion_data = entry_data(mtnme_1, mtnme_2, direct_0, mprester_key)
	entries = entry_ion_data["entries"]
	ion_dict_1 = entry_ion_data["ion_dict_1"]
	if not mtnme_1==mtnme_2:
		ion_dict_2 = entry_ion_data["ion_dict_2"]
	print ion_dict_1
	############################## 1 Element ########################
	if mtnme_1 == mtnme_2:
		ref_state_1=str(ion_dict_1[0]['Reference Solid'])
		ref_dict_1 = {ref_state_1: ion_dict_1[0]['Reference solid energy']}
		entries_aqcorr = aq_correction(entries)

		# #TEMP
		# for i in entries_aqcorr:
		# 	i.correction=0

		stable_solids_minus_h2o = stable_entr(entries_aqcorr)

		pbx_solid_entries = form_e(stable_solids_minus_h2o,
		entries_aqcorr)

		pbx_ion_entries_1 = mke_pour_ion_entr(mtnme_1,
		ion_dict_1, stable_solids_minus_h2o, ref_state_1,
		entries_aqcorr, ref_dict_1)

		all_entries = pbx_solid_entries + pbx_ion_entries_1

		return all_entries

############################## 2 Elements #######################
	else:
		ref_state_1=str(ion_dict_1[0]['Reference Solid'])
		ref_state_2=str(ion_dict_2[0]['Reference Solid'])
		ref_dict_1 = {ref_state_1: ion_dict_1[0]['Reference solid energy']}
		ref_dict_2 = {ref_state_2: ion_dict_2[0]['Reference solid energy']}
		entries_aqcorr = aq_correction(entries)

		# # TEMP
		# for i in entries_aqcorr:
		# 	i.correction=0

		stable_solids_minus_h2o = stable_entr(entries_aqcorr)

		pbx_solid_entries = form_e(stable_solids_minus_h2o,
		entries_aqcorr)

		pbx_ion_entries_1 = mke_pour_ion_entr(mtnme_1,
		ion_dict_1, stable_solids_minus_h2o, ref_state_1,
		entries_aqcorr, ref_dict_1)

		pbx_ion_entries_2 = mke_pour_ion_entr(mtnme_2,
		ion_dict_2, stable_solids_minus_h2o, ref_state_2,
		entries_aqcorr, ref_dict_2)

		all_entries = pbx_solid_entries + pbx_ion_entries_1 + pbx_ion_entries_2

		return all_entries
