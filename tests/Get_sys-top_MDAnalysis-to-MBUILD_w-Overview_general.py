#!/usr/bin/env python

# In[1]:

import sys
import numpy as np
import MDAnalysis as mda
import mbuild as mb
# add the appropriate directory to the path
# get the indir from the current dir
import os
cdir = os.getcwd()
flag = False
wdir = 'systop_utils'
dir_len = len(wdir)
ctr = 0
while (not flag):
    if ( cdir[ctr:dir_len+ctr] == wdir ):
        flag = True
        dir_ind = ctr
    ctr += 1
indir = cdir[:ctr+dir_len-(len(wdir)+1)]
# import inhouse functions
sys.path.append(indir+'/systop_utils/systop_utils/')
import systop_utils

# In[2]:

sys_nm = sys.argv[1]
top_fnm = sys.argv[2]
conf_fnm = sys.argv[3]
top_format = sys.argv[4]
if top_format == 'None' or top_format == 'none' or top_format == 'NONE':
    top_format = None
conf_format = sys.argv[5]
if conf_format == 'None' or conf_format == 'none' or conf_format == 'NONE':
    conf_format = None
universe = mda.Universe(top_fnm, conf_fnm, topology_format=top_format, format=conf_format)

# In[3]:

# Construct System Topology
# mol / fragments
# segments
# residues / monomers
# atoms
system = systop_utils.Compound(name=sys_nm) # store the topology in the mbuild class "Compound"


# In[4]:

# first get the atom attributes from the MDAnalysis universe
# nb - moltypes, molnums, elements, and names are defined according to the information available
# moltypes <--> fragtypes
# molnums <--> fragnums
# elements <--> types
# names <--> types
print('\n')
print('Getting the atom attributes from universe...')
system.get_atom_attributes(universe)

# check that all the expected atom attributes exist
print('system.atoms is '+str(hasattr(system, 'atoms')))
print('system.atoms.moltypes is '+str(hasattr(system.atoms, 'moltypes')))
print(system.atoms.moltypes)
print('system.atoms.molnums is '+str(hasattr(system.atoms, 'molnums')))
print(system.atoms.molnums)
print('system.atoms.resnames is '+str(hasattr(system.atoms, 'resnames')))
print(system.atoms.resnames)
print('system.atoms.resids is '+str(hasattr(system.atoms, 'resids')))
print(system.atoms.resids)
print('system.atoms.elements is '+str(hasattr(system.atoms, 'elements')))
print(system.atoms.elements)
print('system.atoms.names is '+str(hasattr(system.atoms, 'names')))
print(system.atoms.names)


# In[5]:

# use the atom attributes to construct the basic hierarchical structure molecules-->segments-->residues-->atoms
print('\n')
print('Constructing topology...')
system.construct_topology()


# In[6]:

# add the bonds to the topology structure using the MDAnalysis universe
# nb - the mbuild Compound class stores the bonds as a graph (similar to a networkx object)
#    - it is quite expensive to add the bonds one by one...not even sure what the utility of this is, we might just add
#       a reference to the bond list in the force field (or vica versa).
print('Adding bonds from universe...')
system.add_bonds_from_universe(universe)


# In[7]:

# we can add a representative structure to the topology structure
# nb - this could just be a reference to the input/starting configuration
print('Adding positions from universe...')
system.add_positions_from_universe(universe)


# In[9]:


# with the mbuild Compound class, we can easily view the system or subsystems with nglview
#system.children[0].children[0]._visualize_nglview()


# In[10]:


# Here we use the atom attributes to gather a summary/overview of the topology
# nb - this could probably also be done directly from the topology data structure...that would probably be more general/robust
print('\n')
print('Getting the topology summary...')
system.get_topology_summary()


# In[11]:

print('\n')
print('SYSTEM SUMMARY')
print('--------------')

print('The system is '+str(system))


for i_mol, moltype in enumerate(system.molnames):
    print('\n')
    print('MOLTYPE = '+moltype+':'+'\t'+str(system[moltype][0]))
    print('Number of molecules of this type: '+str(system.moltypes[moltype].count))
    if system.moltypes[moltype].n_res == 1:
        print('A dictionary with the counts of each atom element/type for this molecule type: '+str(system.moltypes[moltype].atom_count))
        print('The formula of this molecule type, with respect to the atoms: '+str(system.moltypes[moltype].atom_formula))
        print('\t ATOMS:')
        for atom in system[moltype][0].children:
            print('\t\t'+str(atom))
    else:
        print('A dictionary with the counts of each residue type within this molecule type: '+str(system.moltypes[moltype].residue_count))
        print('The formula of this molecule type, with respect to residues: '+str(system.moltypes[moltype].residue_formula))
        print('The full sequence of residues for this molecule type: '+str(system.moltypes[moltype].residue_sequence))
        for i_res, res in enumerate(system.moltypes[moltype].resnames):
            print('\t RESTYPE = '+res+':'+'\t'+str(system[moltype][0][res][0])) # this is using the first instance of the residue as representative, could be an issue with protein capping groups!
            print('\t A dictionary with the counts of each atom element/type for this residue type (and molecule type): '+str(system.moltypes[moltype].restypes[res].atom_count))
            print('\t The formula of each residue type (for this molecule type), with respect to the atoms: '+str(system.moltypes[moltype].restypes[res].atom_formula))
            print('\t\t ATOMS:')
            for atom in system[moltype][0][res][0].children:
                print('\t\t'+str(atom))