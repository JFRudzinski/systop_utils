import numpy as np
import mbuild as mb

class atoms(object):

    def __init__(self):
        # These are atom attributes -- in practice might just be references to another part of the metadata
        self.moltypes = [] 
        self.molnums = [] 
        self.resnames = [] 
        self.resids = []
        self.elements = [] 
        self.names = []

class restypes(object):

    def __init__(self):
        self.atom_count = None
        self.atom_formula = None

class moltypes(object):

    def __init__(self):
        self.count = None
        self.n_res = None
        self.resnames = None
        self.restypes = {}
        
    def _add_atom_info(self,atom_count=None,atom_formula=None):
        setattr(self, 'atom_count', atom_count)
        setattr(self, 'atom_formula', atom_formula)
        return

    def _add_res_info(self,residue_count=None,residue_formula=None,residue_sequence=None):
        setattr(self, 'residue_count', residue_count)
        setattr(self, 'residue_formula', residue_formula)
        setattr(self, 'residue_sequence', residue_sequence) # not sure if we even need this, only relevant for linear chains
        return 



class Compound(mb.Compound):

    def __init__(self):
        mb.Compound.__init__(self)
        # These are atom attributes -- in practice might just be references to another part of the metadata
        self.atoms = atoms()
        # These summarize the hierarchical structure of the topology
        # atoms, moltypes, restypes could/should be objects
        self.molnames = None
        self.moltypes = {}
        #moltypes()

    def get_atom_attributes(self, universe):
        # set the conventions for the topology based on the available info
        self.atoms.moltypes = get_moltypes(universe)
        for moltype in self.atoms.moltypes:
            self.moltypes[moltype] = moltypes()
        self.atoms.molnums = get_molnums(universe)
        self.atoms.resnames = get_resnames(universe)
        self.atoms.resids = get_resids(universe)
        self.atoms.elements = get_elements(universe)
        self.atoms.names = get_names(universe)

        return

    def construct_topology(self):

        # Generate the hierarchical data structure for the topology
        # molecules --> segments --> residues --> particles
        # nb - MUST call get_atom_attributes() first!
        # nb - Need to generalize for fragments, not just residues
        for i_mol, mol_id in enumerate(np.unique(self.atoms.molnums)): # loop over molecules

            # get the indices of the current molecule
            mol_ids = np.where(self.atoms.molnums == mol_id)[0]
            # get the moltype from the first particle in this molecule
            moltype = self.atoms.moltypes[mol_ids[0]] 
            # add this molecule to the system
            self.add(mb.Compound(),label=moltype+'[$]')
            # get the set of unique resids for this molecule
            mol_resids = np.unique(self.atoms.resids[mol_ids])
            n_res = mol_resids.shape[0] 

            if n_res > 1: # add the residues (and then the particles to the residues)
                self._add_residues(self.children[i_mol], mol_resids, mol_ids)
            elif n_res == 1: # add the particles directly to the molecule
                self._add_particles(self.children[i_mol], mol_ids)

        return 

    def _add_residues(self, compound, mol_resids, mol_ids):
        for i_res, res_id in enumerate(mol_resids): # loop over residues in this molecule
            # get the (global) indices of the current residue
            # nb - constrained within the current molecule (this avoids issues in the case that all atoms are denoted as the same residue, which sometimes happens by default)
            res_ids = np.where(self.atoms.resids[mol_ids] == res_id)[0]
            # get the resname from the first particle in this residue
            resname = self.atoms.resnames[res_ids[0]]
            # add this residue to the molecule
            compound.add(mb.Compound(),label=resname+'[$]')

            self._add_particles(compound.children[i_res], res_ids)

    def _add_particles(self, compound, compound_ids):
        for particle_id in compound_ids: # loop over the particles in this compound
            # create a particle object 
            particle = mb.Particle(name=self.atoms.elements[particle_id])
            # add the object to the residue
            compound.add(particle, label=self.atoms.names[particle_id]+'[$]')

        return


    def add_bonds_from_universe(self, universe):
        # add the bonds to the topology structure
        # This is best of 3 methods I tried, but still rather slow, ~2.5m per 1k bonds on my laptop - faster way?!

        # restrict to 2k bonds for demonstration
        n_bonds_restrict = 2000
        if hasattr(universe, 'bonds'):
            if hasattr(universe.bonds, '_bix'):
                _ = [ self.add_bond((self[int(x[0])],self[int(x[1])])) for x in universe.bonds._bix[:n_bonds_restrict] ]
            else:
                raise ValueError('No attribute "_bix" in universe.bonds')
        else:
            raise ValueError('No attribute "bonds" in the universe.') 
            
        return

    def add_positions_from_universe(self, universe):
        # add a representative configuration from the universe
        if hasattr(universe.atoms, 'positions'):
            self.xyz = universe.atoms.positions
            
        return

    def get_topology_summary(self):

        self.molnames = np.unique(self.atoms.moltypes)

        for moltype in self.molnames: # loop over molecule types

            # MOLTYPE COUNTS
            # get all atom indices associated with this molecule type
            mol_inds = np.where(self.atoms.moltypes == moltype)[0]
            # get the molecule identifier for each atom (i.e., which molecule it belongs to)
            mol_nums = self.atoms.molnums[mol_inds]
            # get the total number of molecules of this type
            self.moltypes[moltype].count = np.unique(mol_nums).shape[0]

            # MOLTYPE RESIDUE FORMULA, COUNTS, and SEQUENCE
            # get the unique set of resids for this molecule
            mol_resids = np.unique(self.atoms.resids[mol_inds])
            # get the (relative) indices associated with the first instance of this molecule type (**indices are relative to the list for this molecule type)
            mol_inds_first_rel = np.where(mol_nums == mol_nums[0])[0]
            # get the (global) indices for the first instance of this molecule type
            mol_inds_first = mol_inds[mol_inds_first_rel]
            # get the corresponding elements (from above)
            elements_mol_first = self.atoms.elements[mol_inds_first]
            # get the number of residues in this molecule
            n_res = np.unique(self.atoms.resids[mol_inds_first]).shape[0]
            # store the number of residues for this molecule
            self.moltypes[moltype].n_res = n_res
            if n_res == 1: # skip the residue level
                # self.moltypes[moltype].residue_count = {}
                # self.moltypes[moltype].residue_formula = {}
                # self.moltypes[moltype].residue_sequence = {}
                # add the summary attributes of the atoms directly to moltypes
                self.moltypes[moltype]._add_atom_info()
                self.moltypes[moltype].atom_count = {}
                # self.moltypes_restypes_atom_formula[moltype] = {}
                self.moltypes[moltype].atom_formula = get_compound_formula(self.moltypes[moltype].atom_count, elements_mol_first)
                continue

            # add the summary attributes of the residues to moltypes
            self.moltypes[moltype]._add_res_info()
            # get the corresponding residue names
            resnames_mol_first = self.atoms.resnames[mol_inds_first]
            # get the corresponding residue identifiers
            resids_mol_first = self.atoms.resids[mol_inds_first]
            # filter for the first instance of each residue, as to not overcount
            # get the counts - number of atoms per residue
            resids_mol_first_uniq, resids_mol_first_cnt = np.unique(resids_mol_first, return_counts=True)
            # get the index of the first atom of each residue
            resids_mol_first_atom_first = np.cumsum(resids_mol_first_cnt)[:-1]
            # add the 0th index manually
            resids_mol_first_atom_first = np.insert(resids_mol_first_atom_first,0,0)
            resnames_mol_first_atom_first = resnames_mol_first[resids_mol_first_atom_first]
            # get the list of unique residues in this molecules, the counts for each, and also the sequence
            self.moltypes[moltype].residue_count = {}
            self.moltypes[moltype].residue_formula = get_compound_formula(self.moltypes[moltype].residue_count, resnames_mol_first_atom_first)
            self.moltypes[moltype].residue_sequence = get_compound_sequence(resnames_mol_first_atom_first)
            

            # This is for getting the atom list for each unique residue within a molecule type
            # get the list of unique residue names for each molecule type
            res_uniq = np.unique(resnames_mol_first)
            # store the unique resnames
            self.moltypes[moltype].resnames = res_uniq
            # self.moltypes[moltype].restypes[res].atom_count = {}
            # self.moltypes[moltype].restypes[res].atom_formula = {}
            for i_res, res in enumerate(res_uniq):
                # add the restype 
                self.moltypes[moltype].restypes[res] = restypes()
                # get the corresponding atom indices
                atom_inds = np.where(resnames_mol_first==res)[0]
                # select only the first instance of this residue
                rel_inds_first_res = np.where(resids_mol_first[atom_inds]==resids_mol_first[atom_inds][0])
                # filter the atom indices
                atom_inds = atom_inds[rel_inds_first_res]
                # get the list of unique atoms in this res and also the counts for each
                self.moltypes[moltype].restypes[res].atom_count = {}
                self.moltypes[moltype].restypes[res].atom_formula = get_compound_formula(self.moltypes[moltype].restypes[res].atom_count, elements_mol_first[atom_inds]) 

        return


def get_fragtypes(universe): 
    # Function for determining the fragment types by hand
    # input - MDAnalysis universe
    # output - array of fragment types, size of n_atoms
    atoms_fragtypes = np.empty(universe.atoms.types.shape, dtype="str")
    ctr_fragtype = 0
    # set the first fragment to type 0
    atoms_fragtypes[universe.atoms.fragments[0]._ix] = ctr_fragtype
    # store the atom types for the first fragment
    frag_unique_atomtypes = []
    frag_unique_atomtypes.append( universe.atoms.types[universe.atoms.fragments[0]._ix] )
    ctr_fragtype += 1
    for i_frag in range(1,universe.atoms.n_fragments): # loop over the remaining fragments
        # current fragment id
        types_i_frag = universe.atoms.types[universe.atoms.fragments[i_frag]._ix]
        # keep track if we find a match for this fragment    
        flag_fragtype_exists = False 
        for j_frag in range(len(frag_unique_atomtypes)-1,-1,-1): # check back through the unique frag types
            types_j_frag = frag_unique_atomtypes[j_frag]
            if len(types_i_frag) != len(types_j_frag): # frags have different number of atoms
                continue
            elif np.all(types_i_frag==types_j_frag): # frags have exactly same list of atom types
                atoms_fragtypes[universe.atoms.fragments[i_frag]._ix] = j_frag
                flag_fragtype_exists = True      

        if not flag_fragtype_exists:
            atoms_fragtypes[universe.atoms.fragments[i_frag]._ix] = ctr_fragtype
            frag_unique_atomtypes.append( universe.atoms.types[universe.atoms.fragments[i_frag]._ix] )
            ctr_fragtype += 1
        
    return atoms_fragtypes

def get_compound_formula(compound_count, children_names):
    # get the list of unique atoms in this res and also the counts for each
    children_count_tup = np.unique(children_names, return_counts=True)
    # concatenate to construct a molecular formula
    children_list = ''
    for child in range(len(children_count_tup[0])):
        # if children_count_tup[1][el] == 1:
        #     children_list += str(children_count_tup[0][child])
        # else:
        children_list += str(children_count_tup[0][child])+'('+str(children_count_tup[1][child])+')'
        compound_count[children_count_tup[0][child]] = children_count_tup[1][child]

    return children_list

def get_compound_sequence(children_names):
    # simply concatenate the list of children
    children_list = children_names[0]
    for child in children_names[1:]:
        children_list += '.'+child

    return children_list

def get_moltypes(universe):
    # assign moltypes between moltypes (gromacs) and fragtypes (other)
    if hasattr(universe.atoms, 'moltypes'): # this is gromacs specific
        return universe.atoms.moltypes
    elif hasattr(universe.atoms, 'fragments'): # otherwise we use fragments to determine the molecules
        atoms_fragtypes = get_fragtypes(universe)
        return atoms_fragtypes
    else:
        raise ValueError('No replacement for moltypes.') 

def get_molnums(universe):
    # assign molnums between molnums (gromacs) and fragindices (other)
    if hasattr(universe.atoms, 'molnums'): # this is gromacs specific
        return universe.atoms.molnums
    elif hasattr(universe.atoms, 'fragindices'): # otherwise we use fragments to determine the molecules
        return universe.atoms.fragindices
    else:
        raise ValueError('No replacement for molnums.') 

def get_resnames(universe):
    if hasattr(universe.atoms, 'resnames'): 
        return universe.atoms.resnames
    elif hasattr(universe.atoms, 'resids'): # otherwise we use resids
        return universe.atoms.resids.astype(str)
    else:
        raise ValueError('No replacement for resnames.')

def get_resids(universe):
    if hasattr(universe.atoms, 'resids'): 
        return universe.atoms.resids
    else:
        raise ValueError('No replacement for resids.')

def get_elements(universe):
    if hasattr(universe.atoms, 'elements'): 
        return universe.atoms.elements
    elif hasattr(universe.atoms, 'types'): # otherwise we use atoms types to name the atoms
        return universe.atoms.types
    else:
        raise ValueError('No replacement for elements.') 

def get_names(universe):
    if hasattr(universe.atoms, 'names'): 
        return universe.atoms.names
    elif hasattr(universe.atoms, 'types'): # otherwise we use atoms types to name the atoms
        return universe.atoms.types
    else:
        raise ValueError('No replacement for names.') 

  
