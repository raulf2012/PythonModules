#!/usr/bin/env python

"""
Author(s): Colin Dickens

TODO
    continue to test
    add reaction printing to __repr__
    add potential/pressure/temperature dependence to __repr__
    add optional keyword arguments to command line interface for setting potential/pressure/temperature
    add Sr, Ru references
"""

#| - Import Modules
import os

import glob
import filecmp

import numpy as np

from ase.atoms import Atoms
from ase.io import read
#__|

class Get_G:
    """
    Class that automatically calculates standard change in free energy between two states by referencing any atoms
    missing between them to a number of gas-phase or aqueous references
    """

    #| - Get_G
    def __init__(self,
        slab,
        ads,
        default_vib_bool=True,
        get_E=False,
        index=-1,
        quiet=False,
        ):
        """
        Get_G(ads,slab) where ads/slab are either paths to traj files or paths to directory containing traj files that have energy and forces.
        The directories that contain these traj files must also contain a calculation directory with pw.inp.
        """
        #| - __init__
        self.default_vib_bool = default_vib_bool
        self.get_E = get_E
        self.quiet = quiet

        if type(slab) == str:
            self.ads_atoms = self.read_atoms(ads,index=index)
            self.slab_atoms = self.read_atoms(slab,index=index)

        elif isinstance(slab, Atoms):

            self.slab_atoms = slab
            self.ads_atoms = ads


        # RF | 181106
        # self.update_params(self.ads_atoms)
        # self.update_params(self.slab_atoms)
        #
        # if self.ads_atoms.PARAMS != self.slab_atoms.PARAMS:
        #
        #     print("%s:"%slab)
        #     print(self.slab_atoms.PARAMS)
        #     print("%s:"%ads)
        #     print(self.ads_atoms.PARAMS)
        #
        #     # print "%s:"%slab
        #     # print self.slab_atoms.PARAMS
        #     # print "%s:"%ads
        #     # print self.ads_atoms.PARAMS
        #
        #     for param in ('dw','pp','pw','xc'):
        #         if self.ads_atoms.PARAMS[param] != self.slab_atoms.PARAMS[param]:
        #             raise Exception("Calculations performed with different parameters")
        #     else:
        #         print("WARNING, ONE CALCULATION IS SPIN-POLARIZED")

        self.update_delta_atoms()
        self.vib_correction()
        self.set_refs()

        # RF | 181106
        # if self.slab_atoms.PARAMS['sp']:
        #     self.compare_magmoms()

        # self.calc_G()
        # self.compare_wf()
        #__|

    def read_atoms(self,path,**kwargs):
        """
        Loads traj file. Path must be direct path to traj files or paths to directory containing them named as either qn.traj or qnXX.traj.
        """
        #| - read_atoms
        if path.find('traj') != -1: #Checks if directory or filename has been specified
            atoms = read(path)
            if '/' in path:
               atoms.PATH = '/'.join(path.split('/')[:-1])
            else:
               atoms.PATH = '.'
        else:
            files = glob.glob(path + '/qn*.traj')
            qn = -1
            for file in files:
                file = file.split('/')[-1]
                if len(file.split(".")[0]) == 2: #qn.traj
                    atoms_path = path + '/qn.traj'
                elif int(file.split('.')[0][2:]) > qn: #qnXX.traj
                    qn = int(file.split('.')[0][2:])
                    atoms_path = path + '/qn%i.traj'%qn
            try:
               atoms = read(atoms_path, **kwargs)
               atoms.PATH = path
            except NameError:
                raise IOError("Could not find traj file associate with " + path)
        if self.fmax(atoms) > 0.05:
            print("WARNING: fmax = %.2f for atoms in %s"%(self.fmax(atoms),atoms.PATH))
        return atoms
        #__|

    def update_params(self,atoms):
        """
        Takes atoms object containing PATH attribute and adds PARAMS dict as attribute with keys pw, dw, xc, pp.
        Assumes PATH contains outdir or calcdir with pw.inp
        """
        #| - update_params
        if os.path.isdir(atoms.PATH + '/outdir'):
            calcdir = atoms.PATH + '/outdir'
        elif os.path.isdir(atoms.PATH + '/calcdir'):
            calcdir = atoms.PATH + '/calcdir'
        else:
            raise IOError('Cannot find calculation directory (outdir or calcdir) for ' + atoms.PATH)

        file = open(calcdir + '/pw.inp')
        lines = file.readlines()
        file.close()

        self.params = {}
        self.params['sp'] = False

        for line in lines:
            if line[2:9] == 'ecutwfc':
                self.params['pw'] = int(float(line.split('=')[-1][:-4])*rydberg)
            if line[2:9] == 'ecutrho':
                self.params['dw'] = int(float(line.split('=')[-1][:-4])*rydberg)
            if line[2:11] == 'input_dft':
                self.params['xc'] = line.split('=')[-1][1:-3]
            if line[2:12] == 'pseudo_dir':
                self.params['pp'] = line.split('=')[-1][1:-3]
            if line[2:8] == 'nspin=':
                self.params['sp'] = True


        for key in synonyms:
           if self.params[key] in synonyms[key]:
               self.params[key] = synonyms[key][self.params[key]]

        atoms.PARAMS = self.params
        #__|

    def update_delta_atoms(self):
        """
        Update dictionary self.data_atoms with difference in chemical formulas between ads and slab.
        Positive numbers represent species that ads has and slab doesn't.
        """
        #| - update_delta_atoms
        self.delta_atoms = {}
        ads_syms = self.ads_atoms.get_chemical_symbols()
        slab_syms = self.slab_atoms.get_chemical_symbols()
        ads_dict = {}
        slab_dict = {}

        for sym in ads_syms:
            if sym in ads_dict:
                ads_dict[sym] += 1
            else:
                ads_dict[sym] = 1

        for sym in slab_syms:
            if sym in slab_dict:
                slab_dict[sym] += 1
            else:
                slab_dict[sym] = 1

        for sym in ads_dict:
            try:
                self.delta_atoms[sym] = ads_dict[sym] - slab_dict[sym]
            except KeyError: #no sym in slab_dict
                self.delta_atoms[sym] = ads_dict[sym]

        for sym in slab_dict:
            if sym not in self.delta_atoms:
                self.delta_atoms[sym] = -slab_dict[sym]

        print(self.delta_atoms)
        print("LKSDJFKLDJS_------_-_-__----")

        # for key in self.delta_atoms.keys():
        for key in list(self.delta_atoms):
            if self.delta_atoms[key] == 0:
                del self.delta_atoms[key]
        #__|

    def vib_correction(self):
        """
        Attempt to add explicitly calculated vibrational correction from vib directory within PATH of ads_atoms and slab_atoms.
        Otherwise, default vibrational corrections are used by calling default_vib()
        """
        #| - vib_correction

        def vib_indices(atoms):
            """
            Return a dict with symbol/count key/value pairs given an atoms object with PATH attribute
            """
            #| - vib_indices
            dict = {}

            vib_pkls = glob.glob(atoms.PATH + '/vib/vib.*.pckl')
            indices = []
            for vib_pkl in vib_pkls:
                if vib_pkl != atoms.PATH + '/vib/vib.eq.pckl':
                    i = vib_pkl.split('/')[-1].split('.')[1][:-2]
                    if int(i) not in indices:
                        indices.append(int(i))

            for index in indices:
                sym = atoms[index].symbol
                if sym in dict:
                    dict[sym] +=1
                else:
                    dict[sym] = 1

            return dict
            #__|

        def parse_corr(path):
            """
            Return vibrational correction from myjob.out in path. Return -1 if myjob.out cannot be read (imag freqs?)
            """
            #| - parse_corr
            file = open(path)
            lines = file.readlines()
            file.close()
            for line in lines:
                if len(line.split()) > 1:
                    if line.split()[0] == 'G' or line.split()[0] == 'F':
                        corr = float(line.split()[1])
                        return corr
            return -1

            """
            try:
                try:
                   corr = float([item for item in [line for line in lines if 'G' in line][0].split(' ') if item.strip()][1])
                   return corr
                except:
                   corr = float([item for item in [line for line in lines if 'F' in line][0].split(' ') if item.strip()][1])
                   return corr
            except:
                return -1
            """
            #__|

        if self.get_E:
            self.ads_corr = 0
            self.slab_corr = 0
            return

        if self.default_vib_bool:
            self.default_vib()
            return

        ads_dict = vib_indices(self.ads_atoms)
        slab_dict = vib_indices(self.slab_atoms)

        self.delta_vib = {}

        for sym in ads_dict:
            try:
                self.delta_vib[sym] = ads_dict[sym] - slab_dict[sym]
            except KeyError: #no sym in slab_dict
                self.delta_vib[sym] = ads_dict[sym]

        for sym in slab_dict:
            if sym not in self.delta_vib:
                self.delta_vib[sym] = -slab_dict[sym]

        for sym in self.delta_vib.keys():
            if self.delta_vib[sym] == 0:
                del self.delta_vib[sym]


        if self.delta_vib == self.delta_atoms: #explicit vibrational corrections appear valid
            if ads_dict != {}:
                self.ads_corr = parse_corr(self.ads_atoms.PATH + '/vib/myjob.out')
                if self.ads_corr == -1:
                    print("Correct number of vibrational modes, but cannot read myjob.out (imag freqs?)")
                    self.default_vib()
                    return
            else:
                self.ads_corr = 0

            if slab_dict != {}:
                self.slab_corr = parse_corr(self.slab_atoms.PATH + '/vib/myjob.out')
                if self.slab_corr == -1:
                    print("Correct number of vibrational modes, but cannot read myjob.out (imag freqs?)")
                    self.default_vib()
                    return
            else:
                self.slab_corr = 0
        else:
            self.default_vib()
            return
        #__|

    def default_vib(self):
        """
        Calculate vibrational corrections using estimates based on atom identity.
        """
        #| - default_vib
        self.default_vib_bool = True
        self.ads_corr = 0
        self.slab_corr = 0
        for sym in self.delta_atoms:
            if self.delta_atoms[sym] > 0:
                self.ads_corr += default_vib_dict[sym]*self.delta_atoms[sym]
            else:
                self.slab_corr -= default_vib_dict[sym]*self.delta_atoms[sym]
        #__|

    def set_refs(self):
        """
        Formulate references for atoms in self.delta_atoms. Currently chooses reference indicated as default.
        """
        #| - set_refs
        self.references = {}
        for sym in self.delta_atoms:
            for atom in references:
                if atom == sym:
                    for ref in references[atom]:
                        if len(references[atom][ref]) == 3: #use default reference
                            self.references[sym] = references[atom][ref][:2]
        #__|

    def calc_G(self):
        """
        Perform free energy calculation at standard conditions and record thermodynamic dependencies on pressure, tempreature, potential, etc.
        """
        #| - calc_G
        self.G_std = self.ads_atoms.get_potential_energy() - self.slab_atoms.get_potential_energy()
        self.G_std += self.ads_corr - self.slab_corr
        self.G_thermo = {}
        for sym in self.delta_atoms:
            for DFT_ref in self.references[sym][0]:
                n = DFT_ref[1]
                DFT_G = self.get_DFT_ref(DFT_ref[0])
                self.G_std -= self.delta_atoms[sym]*DFT_G*n
            for thermo,n in self.references[sym][1]:
                if thermo in self.G_thermo:
                    self.G_thermo[thermo] -= n*self.delta_atoms[sym]
                else:
                    self.G_thermo[thermo] = -n*self.delta_atoms[sym]
        #__|

    def get_DFT_ref(self,ref):
        """
        Pull appropriate DFT reference energy from database according to computational parameters
        """
        #| - get_DFT_ref
        xc = self.params['xc']
        pwdw =(self.params['pw'],self.params['dw'])
        pp = self.params['pp']

        if self.get_E:
            DFT_references = DFT_E_references
        else:
            DFT_references = DFT_G_references

        if xc in DFT_references[ref]:
            if pwdw in DFT_references[ref][xc]:
                if pp in DFT_references[ref][xc][pwdw]:
                    return DFT_references[ref][xc][pwdw][pp]
                else:
                    for pp_ref in DFT_references[ref][xc][pwdw]:
                        if self.compare_pp(pp,pp_ref,DFT_references[ref]['syms']):
                            return DFT_references[ref][xc][pwdw][pp_ref]

        raise Exception("No reference found for %s with %s @ %s with %s"%(ref,xc,pwdw,pp))
        #__|

    def compare_pp(self,pp1,pp2,syms):
        """
        """
        #| - compare_pp
        for sym in syms:
            if not filecmp.cmp("%s/%s.UPF"%(pp1,sym),"%s/%s.UPF"%(pp2,sym)):
                return False
        return True
        #__|

    def fmax(self,atoms):
        """
        """
        #| - fmax
        forces = atoms.get_forces()
        max = 0
        for force in forces:
             tot = (force[0]**2 + force[1]**2 + force[2]**2)**0.5
             if tot > max:
                 max = tot
        return max
        #__|

    def compare_magmoms(self):
        """
        """
        #| - compare_magmoms
        def nearest_atom(atoms,position):
            "Returns atom nearest to position"
            position = np.array(position)
            dist_list = []
            for atom in atoms:
                dist = np.linalg.norm(position - atom.position)
                dist_list.append(dist)

            return atoms[np.argmin(dist_list)]

        if len(self.ads_atoms) >= len(self.slab_atoms):
            ads = self.ads_atoms
            slab = self.slab_atoms
            indexed_by = "slab"
            not_indexed_by = "ads"
        else:
            slab = self.ads_atoms
            ads = self.slab_atoms
            indexed_by = "ads"
            not_indexed_by = "slab"

        delta_magmoms = []
        ads_indices_used = []
        for atom in slab:
            ads_atom = nearest_atom(ads,atom.position)
            if not self.quiet:
               if ads_atom.symbol != atom.symbol: print("WARNING! MAGMOM COMPARISON FAILURE")
            ads_indices_used.append(ads_atom.index)
            delta_magmoms.append(atom.magmom - ads_atom.magmom)

        ads_indices_not_used = []
        for i in range(len(ads)):
            if i not in ads_indices_used:
                ads_indices_not_used.append(i)

        # RF | 181106
        # self.delta_magmoms = zip(range(len(slab)), delta_magmoms)
        self.delta_magmoms = list(zip(range(len(slab)), delta_magmoms))
        self.delta_magmoms.sort(key=lambda x: abs(x[1]),reverse=True)

        common = ""
        uncommon = ""
        for i in range(8):
            atom = slab[self.delta_magmoms[i][0]]
            common += "%s%d: %.2f\t"%(atom.symbol,atom.index,self.delta_magmoms[i][1])
        for i in ads_indices_not_used:
            uncommon += "%s%d: %.2f\t"%(ads[i].symbol,ads[i].index,ads[i].magmom)

        if self.quiet:
            return
        else:

            print("~"*6 + "MAGNETIC MOMENT COMPARISON" + "~"*6)
            print("Largest magnetic moment discrepancies (indexed by %s)"%indexed_by)
            print(common)
            print("Magnetic moments only present in %s"%not_indexed_by)
            print(uncommon + "\n")

            # print("~"*6 + "MAGNETIC MOMENT COMPARISON" + "~"*6)
            # print("Largest magnetic moment discrepancies (indexed by %s)"%indexed_by)
            # print(common)
            # print("Magnetic moments only present in %s"%not_indexed_by)
            # print(uncommon + "\n")
        #__|

    def rms_displacement(self):
        """
        """
        #| - rms_displacement
        displacements = np.zeros(np.min((len(self.slab_atoms),len(self.ads_atoms))))
        for i in range(len(displacements)):
            displacements[i] = np.linalg.norm(self.slab_atoms[i].position-self.ads_atoms[i].position)**2
        return np.sqrt(displacements.mean())
        #__|

    def compare_wf(self):
        """
        """
        #| - compare_wf
        self.wf_slab = self.read_wf(self.slab_atoms.PATH)
        self.wf_ads = self.read_wf(self.ads_atoms.PATH)

        if 'N/A' in [self.wf_slab,self.wf_ads]:
           self.d_wf = 'N/A'
        else:
           self.d_wf = self.wf_ads - self.wf_slab
           self.avg_wf = (self.wf_ads + self.wf_slab)/2
           self.d_mu = 0.0055 * self.d_wf * np.linalg.norm(np.cross(self.ads_atoms.cell[0],self.ads_atoms.cell[1]))
        #__|

    def read_wf(self,path):
        """
        """
        #| - read_wf
        try:
           f = open(path + '/out.WF')
        except:
            return 'N/A'

        line = f.readlines()[0]
        return float(line.split(',')[0][1:])
        #__|


    #| - __old__ | old __repr__
    # def __repr__(self):
    #     """
    #     String representation of free energy calculation. Example output:
    #
    #     Ir48O100 + H+/e- --> HIr48O100   dG = (-0.5 + eU_RHE)
    #     pw/dw = 600/6000    xc = RPBE   pp = /home/vossj/suncat/esdld/psp
    #     Default vibrational corrections applied to adsorbates
    #     Other possible references include H2...
    #     """
    #     #| - __repr__
    #     string = ""
    #     if self.get_E:
    #         string += "dE = %.3f eV"%self.G_std
    #     else:
    #         string += "dG = %.3f eV"%self.G_std
    #
    #         for thermo in self.G_thermo:
    #             if self.G_thermo[thermo] > 0:
    #                string += " + %d%s"%(self.G_thermo[thermo],thermo)
    #             elif self.G_thermo[thermo] < 0:
    #                string += " - %d%s"%(-self.G_thermo[thermo],thermo)
    #
    #     string += '\n'
    #     string += "xc = %s\tpw/dw = %d/%d\tpp = %s"%(self.params['xc'],self.params['pw'],self.params['dw'],self.params['pp'])
    #
    #     if not self.get_E:
    #         if self.default_vib_bool:
    #             string += "\nUsing default vibrational corrections for adsorbates"
    #         else:
    #             string += "\nUsing vibrational corrections for adsorbates found in directory"
    #
    #     string += "\nRMS Displacement = %.3f Angstroms"%self.rms_displacement()
    #     if type(self.d_wf) == float:
    #        string += "\nIS WF = %.2f eV, FS WF = %.2f eV, Change in WF = %.2f eV, Change in Dipole = %.2f eA, Average WF = %.2f eV"\
    #                %(self.wf_slab,self.wf_ads,self.d_wf,self.d_mu,self.avg_wf)
    #     return string
    #     #__|
    #__|

    #__|







#| - out_of_sight

rydberg = 13.6057 #rydberg to eV conversion

default_vib_dict = {'H':0.3, 'O':0.05, 'N': 0.05, 'Ru':0.}

references = {
        'H':{
            'H2':([('H2',0.5)],[('kB*T*ln(P_H2)',0.5)]),
            'H+/e-':([('H2',0.5)],[('e*U_RHE',-1)],'DEF')
            },
        'O':{
            'O2':([('H2O(l)',1),('H2',-1),(2.46)],[('kB*(300 K)*ln(P_O2)',0.5)]),
            'H2O(l)':([('H2O(l)',1),('H2',-1)],[('e*U_RHE',2)],'DEF')
            },
        'N':{
            'N2':([('N2',0.5)],[('kB*T*ln*(P_N2)',0.5)],'DEF')
            },
        'Ru':{
            'RuO4':([('RuO4(l)',1),('H2O(l)',-4),('H2',4)],[('e*U_RHE',-8)],'DEF')
            },
        'Li':{
            'Li':([('Li',1)],[('N/A',0)],'DEF')
            },
        }


#Preliminary dictionary of DFT reference chemical potentials (assuming ideal gas or electronic energy for solids)
DFT_G_references = {
        'H2':{
            'RPBE':{
                (600,6000):{
                    '/home/vossj/suncat/esdld/psp':-32.1129703764,
                    '/scratch/users/colinfd/psp/gbrv':-32.0745195582
                    },
                (500,5000):{
                    '/home/vossj/suncat/esdld/psp':-32.097274273
                    }
                },
            'PBE':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-31.7495756638,
                    '/home/vossj/suncat/esdld/psp':-31.7878744197
                    }
                },
            'BEEF':{
                (600,6000):{
                    '/home/vossj/suncat/esdld/psp':-33.0011548009,
                    '/scratch/users/colinfd/psp/gbrv':-32.9557311734,
                    },
                (500,5000):{
                    '/home/vossj/suncat/esdld/psp':-32.9825919425,
                    '/scratch/users/colinfd/psp/gbrv':-32.955410852
                    },
                (550,5500):{
                    '/home/vossj/suncat/esdld/psp':-32.9946770046,
                    }
                },
            'syms':['H'] #used for comparing pseudopotentials
            },
        'H2O(l)':{
            'RPBE':{
                (600,6000):{
                    '/home/vossj/suncat/esdld/psp':-491.396166885,
                    '/scratch/users/colinfd/psp/gbrv':-472.461551626
                    },
                (500,5000):{
                    '/home/vossj/suncat/esdld/psp':-491.260712609
                    }
                },
            'PBE':{
                (600,6000):{
                    '/home/vossj/suncat/esdld/psp':-490.525317314,
                    '/scratch/users/colinfd/psp/gbrv':-471.663372251
                    }
                },
            'BEEF':{
                (600,6000):{
                    '/home/vossj/suncat/esdld/psp':-496.410477954,
                    '/scratch/users/colinfd/psp/gbrv':-476.619940391,
                    },
                (500,5000):{
                    '/home/vossj/suncat/esdld/psp':-496.404846831
                    }
                },
            'syms':['H','O'] #used for comparing pseudopoentials
            },
        'N2':{
            'BEEF':{
                (500,5000):{
                    '/scratch/users/colinfd/psp/gbrv':-553.966954842,
                    '/home/vossj/suncat/esdld/psp':-555.426924645
                    }
                },
            'syms':['N'] #used for comparing pseudopoentials
            },
        'RuO4(l)':{
            'RPBE':{
                (600,6000):{
                    '/home/vossj/suncat/esdld/psp':-2512.43779553,
                    '/scratch/users/colinfd/psp/esp_psp_w_gbrvRu':-4475.37701007
                    }
                },
            'BEEF':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/esp_psp_w_gbrvRu':-4503.54664532
                    }
                },
            'syms':['Ru','O']
            }
        }


DFT_E_references = {
        'H2':{
            'BEEF':{
                (500,5000):{
                    '/scratch/users/colinfd/psp/gbrv':-32.9198000182,
                    '/home/vossj/suncat/esdld/psp':-32.9423056024,
                    },
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-32.9199219252,
                    '/home/vossj/suncat/esdld/psp':-32.9629591797,
                    }
                },
            'RPBE':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-32.0316308866,
                    '/home/vossj/suncat/esdld/psp':-32.0667357502
                    },
                },
            'PBE':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-31.7044722137,
                    '/home/vossj/suncat/esdld/psp':-31.7393674136
                    }
                },
            'syms':['H'] #used for comparing pseudopoentials
            },
        'N2':{
            'BEEF':{
                (500,5000):{
                    '/scratch/users/colinfd/psp/gbrv':-553.610425617,
                    '/home/vossj/suncat/esdld/psp':-555.069481444
                    }
                },
            'syms':['N'] #used for comparing pseudopoentials
            },

        'H2O(l)':{
            'RPBE':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-472.468333263,
                    '/home/vossj/suncat/esdld/psp':-491.389771898
                    }
                },
            'PBE':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-471.668649297,
                    '/home/vossj/suncat/esdld/psp':-490.517735036
                    }
                },
            'BEEF':{
                (600,6000):{
                    '/scratch/users/colinfd/psp/gbrv':-476.630765091,
                    '/home/vossj/suncat/esdld/psp':-496.411376542
                    }
                },

            'syms':['O','H'] #used for comparing pseudopotentials
            },

        'Li':{
             'BEEF':{
                 (500,5000):{
                     '/home/vossj/suncat/esdld/psp':-204.787513924
                     },
                 },
             }
        }

##synonyms
synonyms = {
        'pp':{
            #2014 gbrv
            '/nfs/slac/g/suncatfs/colinfd/psp/esp_psp_w_gbrvRu':
            '/scratch/users/colinfd/psp/esp_psp_w_gbrvRu',

            '/nfs/slac/g/suncatfs/colinfd/psp/gbrv':
            '/scratch/users/colinfd/psp/gbrv',

            '/home/vossj/suncat/psp/gbrv1.5pbe':
            '/scratch/users/colinfd/psp/gbrv',

            '/global/project/projectdirs/m2997/colinfd/psp/gbrv':
            '/scratch/users/colinfd/psp/gbrv',

            #v2 (default QE)
            '/nfs/slac/g/suncatfs/sw/external/esp-psp/v2':
            '/home/vossj/suncat/esdld/psp',

            '/nfs/slac/g/suncatfs/sw/rh6/external/../../external/esp-psp/v2':
            '/home/vossj/suncat/esdld/psp',

            '/global/project/projectdirs/m2997/colinfd/psp/esp_psp/':
            '/home/vossj/suncat/esdld/psp',

            '/scratch/users/colinfd/psp/esp_psp':
            '/home/vossj/suncat/esdld/psp',
            },
        'xc':{
           'BEEF-vdW':'BEEF'
        }
}

#__|

# if __name__ == "__main__":
#     import argparse
#     parser = argparse.ArgumentParser()
#     parser.add_argument('slab', type=str)
#     parser.add_argument('ads', type=str)
#     parser.add_argument('-vib', '--vibrations', action = 'store_false', help='Look for vibrational calculations in vib/')
#     parser.add_argument('-E', '--elec_only', action = 'store_true', help='Calculate change in electronic energy (no TS, ZPE, etc)')
#     parser.add_argument('-i', '--index', type=int, help = 'Index to use for both traj files', default = -1)
#
#     args = parser.parse_args()
#     G = Get_G(args.slab,args.ads,default_vib_bool = args.vibrations, get_E = args.elec_only, index=args.index)
#     print G
