# -*- coding: utf-8 -*-

def entry_data(mtnme_1, mtnme_2, direct_0, mprester_key):
    """
    Obtaining entry and ion data from local source. If the data cannot be
    found, it will do an HTTP request to the Materials Project database.

    Args:
        mtnme_1: Material 1
        mtnme_2: Material 2
        direct_0: Directory of entry database for binary systems
    """
    # | -  entry_data
    from pymatgen.matproj.rest import MPRester # Access species in MP
    import warnings; import json; from monty.json import MontyDecoder
    warnings.filterwarnings('ignore')    # Ignore errors related to HTTP request
    mpr = MPRester(mprester_key)    # REST api adaptor, INPUT

    try:
        direct_0 = direct_0        #INPUT#
        direct = direct_0+mtnme_1+'_'+mtnme_2+'/'
        entries = json.loads(open(direct+'mp_entries.txt','r').read(), cls=MontyDecoder)
        ion_dict_1 = json.loads(open(direct+'ion_data_1.txt','r').read())

        if not mtnme_1==mtnme_2:
            ion_dict_2 = json.loads(open(direct+'ion_data_2.txt','r').read())
    except:
        pass

################################################################################
    ## Obtains all entries in MP corresponding to input element, O and H
    if mtnme_1==mtnme_2:
        try:
            entries = entries
        except:
            print 'pd_make.entry_data - made http request line 64'
            entries = mpr.get_entries_in_chemsys([mtnme_1, 'O', 'H'])

        ## Dictionary of reference state:experimental formation energy
        try:
            ion_dict_1 = ion_dict_1
        except:
            print 'pd_make.entry_data - made http request line 69'
            ion_dict_1 = mpr._make_request('/pourbaix_diagram/reference_data/'+mtnme_1)
        out_dict = {"entries" : entries, "ion_dict_1" : ion_dict_1}
        return out_dict

    if not mtnme_1==mtnme_2:
        try:
            entries = entries
        except:
            print 'pd_make.entry_data - made http request - line 146'
            ## Obtains all entries in MP corresponding to input element, O and H
            entries = mpr.get_entries_in_chemsys([mtnme_1, mtnme_2, 'O', 'H'])
            ## Dictionary of reference state:experimental formation energy
        try:
            ion_dict_1 == ion_dict_1
            ion_dict_2 == ion_dict_2
        except:
            print 'pd_make.entry_data - made http request - line 154'
            ion_dict_1 = mpr._make_request('/pourbaix_diagram/reference_data/'+mtnme_1)
            ion_dict_2 = mpr._make_request('/pourbaix_diagram/reference_data/'+mtnme_2)
        out_dict = {"entries" : entries, "ion_dict_1" : ion_dict_1, "ion_dict_2" : ion_dict_2}
        return out_dict
        # __|

def remove_duplicate_entries(entry_list):
    """
    """
    # | -  remove_duplicate_entries
    def contains_entry(entry_list, ent):
        """
        Helpful to filter duplicate entries, if entry
        is in entry_list, return True
        Args:
            entry_list: list of pymatgen entries to consider
            entry: entry which will be analyzed
        """

        ent_id = ent.entry_id
        ent_E = ent.energy_per_atom
        ent_redfor = ent.composition.reduced_formula
        for e in entry_list:
            if e.entry_id == ent_id or (abs(ent_E - e.energy_per_atom) < 1e-6
            and ent_redfor == e.composition.reduced_formula):
                return True

    entry_list_new = list()
    for entry in entry_list:
        if not contains_entry(entry_list_new, entry):
            entry_list_new.append(entry)

    return entry_list_new
    # __|

def aq_correction(entries):
    """
    Applies the Materials Project Aqueous Compatibility scheme for mixing GGA
    and GGA+U to a list of entries.
    Removes entries which aren't compatible with the mixing scheme

    Args:
        entries: List of entries on which the correction will be applied
    """
    # | -  aq_correction
    from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility

    def contains_entry(entry_list, ent):
        """
        Helpful to filter duplicate entries, if entry
        is in entry_list, return True
        Args:
            entry_list: list of pymatgen entries to consider
            entry: entry which will be analyzed
        """

        ent_id = ent.entry_id
        ent_E = ent.energy_per_atom
        ent_redfor = ent.composition.reduced_formula
        for e in entry_list:
            if e.entry_id == ent_id or (abs(ent_E - e.energy_per_atom) < 1e-6
            and ent_redfor == e.composition.reduced_formula):
                return True

    aqcompat = MaterialsProjectAqueousCompatibility() #Implements the GGA/GGA+U mixing scheme,

    entries_aqcorr = list()
    for entry in entries:
        aq_corrected_entry = aqcompat.process_entry(entry) #Corrections, if none applicable, gets rid of entry

        if not contains_entry(entries_aqcorr, aq_corrected_entry): #If entry already in entries_aqcorr don't add
            entries_aqcorr.append(aq_corrected_entry)

    return entries_aqcorr
# __|

def stable_entr(entries_aqcorr):
    """
    Evaluate a entries in list for stability and discard unstable entries
    Remove H2, O2, H2O, and H2O2 from the list
    Calculate the formation using the species in the chemical ensemble

    Args:
        entries_aqcorr: List of entries, usually they have been run through the
        aqueous compatibility module
    """
    # | -  stable_entr
    from pymatgen.analysis.pourbaix.entry import PourbaixEntry
    from pymatgen.phasediagram.maker import PhaseDiagram

    ## Generate phase diagram to consider only solid entries stable in water
    pd = PhaseDiagram(entries_aqcorr)
    stable_solids = pd.stable_entries
    stable_solids_minus_h2o = [entry for entry in stable_solids if
        entry.composition.reduced_formula not in ["H2", "O2", "H2O", "H2O2"]]

    return stable_solids_minus_h2o
    # __|

def form_e(stable_solids_minus_h2o, entries_aqcorr, gas_gibbs=True):
    """
    Calculate the formation energy for the entries in stable_solids_minus_h2o
    Reduce by stoicheometric factor if applicable (ex. Fe4O4)

    Args:
        stable_solids_minus_h2o: list of stable solids without O2, H2O, H2, and
    H2O2 (from stable_entr)
        entries_aqcorr: entries list before being modified by stable_entr
    """
    # | -  form_e

    # | -  Imported Modules
    from pymatgen.analysis.pourbaix.entry import PourbaixEntry
    from pymatgen.phasediagram.maker import PhaseDiagram
    from entry_methods import base_atom
    from energy_scheme import ref_atoms_dict
    # __|

    ref_atom_energies = ref_atoms_dict()
    e_o = ref_atom_energies['e_o']
    e_h = ref_atom_energies['e_h']

    pd = PhaseDiagram(entries_aqcorr)

    base_at = base_atom(stable_solids_minus_h2o)
    d = {}
    for i in stable_solids_minus_h2o:
        if i.name in base_at:
            d[i.name] = i.energy_per_atom

    def form_energy(entry, solid_ref_energy_dict, gas_gibbs=True):
        """
        Calculates the formation energy of an entry relative to the VASP
        energies of the reference atoms. Gibbs/entropy corrections for the gas
        phase are optionable

        Args:
            entry: Entry whose formation energy will be calculated
            solid_ref_energy_dict: Dictionary of reference atom energies
            gas_gibbs: Whether the gas phase reference atoms will have an
            entropic correction
        """
        if gas_gibbs==True:
            ref_dict = {}
            ref_dict['O'] = e_o-0.3167
            ref_dict['H'] = e_h-0.36218
        elif gas_gibbs==False:
            ref_dict = {}
            ref_dict['O'] = e_o
            ref_dict['H'] = e_h

        z = ref_dict.copy()
        z.update(solid_ref_energy_dict)
        ref_dict = z

        elem_dict = entry.composition.get_el_amt_dict()
        entry_e = entry.energy

        for elem in entry.composition.elements:
            elem_num = elem_dict[elem.symbol]
            entry_e = entry_e - elem_num*ref_dict[elem.symbol]
        return entry_e

    pbx_solid_entries = []
    for entry in stable_solids_minus_h2o:
        pbx_entry = PourbaixEntry(entry)
        pbx_entry.g0_replace(form_energy(entry, d, gas_gibbs=gas_gibbs)) #Replace E with form E relative to ref elements
        # pbx_entry.g0_replace(pd.get_form_energy(entry)) #Replace E with form E relative to ref elements
        pbx_entry.reduced_entry() #Reduces parameters by stoich factor (ex. Li2O2 -> LiO)
        pbx_solid_entries.append(pbx_entry)

    return pbx_solid_entries
    # __|

###############################################################################
#                    ██  ██████  ███    ██ ███████
#                    ██ ██    ██ ████   ██ ██
#                    ██ ██    ██ ██ ██  ██ ███████
#                    ██ ██    ██ ██  ██ ██      ██
#                    ██  ██████  ██   ████ ███████
###############################################################################

def ref_entry_find(stable_solids_minus_h2o, ref_state):
    """
    Args:
        stable_solids_minus_h2o:
        ref_state:
    """
    # | -  ref_entry_find

    # print ref_state+'_______________'    # TEMP
    #Chk# Ion solid material reference state (defined by the dict 'Reference Solid' entry in ion_dict)
    # is in the generated material list (if YES assign to ref_entry var), if NO generate error
    for entry in stable_solids_minus_h2o:
        # print entry.composition.reduced_formula+'#_#_'    # TEMP
        if entry.composition.reduced_formula==ref_state: ref_entry=entry; break
        else: ref_entry=[]
    #if not ref_entry: return '05 - Error with '+ref_state+' solid reference data for:'+mtnme
    if not ref_entry:
        print '05 - Error with '+ref_state+' solid reference data'
        return '05 - Error with '+ref_state+' solid reference data'

    return ref_entry
    # __|

def ref_entry_stoich(ref_entry):
    """
    """
    # | -  ref_entry_stoich
    ref_stoich_fact=ref_entry.composition.get_reduced_composition_and_factor()[1]
    return ref_stoich_fact
    # __|

def mke_pour_ion_entr(mtnme, ion_dict, stable_solids_minus_h2o, ref_state, entries_aqcorr, ref_dict):
    """
    """
    # | -  mke_pour_ion_entr
    from pymatgen import Element # Accesses properties of element
    from pymatgen.core.ion import Ion
    from pymatgen.phasediagram.maker import PhaseDiagram
    from pymatgen.analysis.pourbaix.entry import PourbaixEntry, IonEntry

    from pd_make import ref_entry_find, ref_entry_stoich
    pd = PhaseDiagram(entries_aqcorr)

    ref_entry = ref_entry_find(stable_solids_minus_h2o, ref_state)
    ref_stoich_fact = ref_entry_stoich(ref_entry)

    ## Calculate DFT reference E for ions (Persson et al, PRB (2012))
    dft_for_e=pd.get_form_energy(ref_entry)/ref_stoich_fact # DFT formation E, normalized by composition "factor"
    ion_correction_1 = dft_for_e-ref_dict[ref_state] # Difference of DFT form E and exp for E of reference

    el = Element(mtnme)
    pbx_ion_entries_1 = []
    for id in ion_dict:
        comp = Ion.from_formula(id['Name'])     # Ion name-> Ion comp name (ex. Fe[3+] -> Ion: Fe1 +3)
                               # comp.composition[el] : number of Fe atoms in ion
        num_el_ref = (ref_entry.composition[el]) / ref_stoich_fact # number of element atoms in reference
        factor = comp.composition[el] / num_el_ref    # Stoicheometric factor for ionic correction
                                   # (i.e. Fe2O3 ref but calc E for Fe[2+] ion)
        energy = id['Energy'] + ion_correction_1 * factor

        #TEMP_PRINT
        if id['Name']=='Pd[2+]':
            energy = 123
        # print id['Name']

        pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy))
        pbx_entry_ion.name = id['Name']
        pbx_ion_entries_1.append(pbx_entry_ion)

    return pbx_ion_entries_1
    # __|
