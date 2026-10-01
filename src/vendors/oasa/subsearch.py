#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------

from . import smiles
from .graph.digraph import digraph

import os

class substructure_search_manager( object):

  substructure_def_file = os.path.join( os.path.dirname( __file__), "subsearch_data.txt")
  ring_def_file = os.path.join( os.path.dirname( __file__), "subsearch_rings.txt")

  def __init__( self):
    self.structures = digraph()  # graph describing the relations between individual structures
    self.search_trees = []  # list of substructure instances that for a tree
    self.rings = {}
    self.fill_data()

  def fill_data( self):
    from . import subsearch_data
    for struct in subsearch_data.structures:
      compound_type, name, smiles_string, to_ignore = struct
      sub = substructure( name, compound_type, smiles=smiles_string, atoms_to_ignore=to_ignore)
      v = self.structures.create_vertex()
      v.value = sub
      self.structures.add_vertex( v)
    for ring_desc in subsearch_data.rings:
      name, smiles_string, ring_hash = ring_desc
      rng = ring( name, smiles_string, ring_hash=ring_hash)
      self.rings[ ring_hash] = rng
    # the following takes some time
    self._analyze_structure_dependencies()
    

  def _read_structure_file( self, name=""):
    """may be used to read data directly from source txt files.
    is deprecated in favor of automatically build subsearch_data.py module"""
    f = file( name or self.substructure_def_file, "r")
    for line in f:
      l = line.strip()
      if l and not l.startswith( "#"):
        parts = l.split(";")
        if len( parts) < 3:
          raise ValueError( "wrong line in data file: '%s'" % l)
        compound_type, name, smiles_string = [x.strip() for x in parts[:3]]
        if len( parts) > 3:
          to_ignore = parts[3].strip()
        else:
          to_ignore = ""
        to_ignore = map( int, filter( None, to_ignore.split(",")))
        if not name.strip():
          name = compound_type
        sub = substructure( name, compound_type, smiles=smiles_string, atoms_to_ignore=to_ignore)
        v = self.structures.create_vertex()
        v.value = sub
        self.structures.add_vertex( v)
    self._analyze_structure_dependencies()

  def _read_ring_file( self, name=""):
    """may be used to read data directly from source txt files.
    is deprecated in favor of automatically build subsearch_data.py module"""
    f = file( name or self.ring_def_file, "r")
    for line in f:
      l = line.strip()
      if l and not l.startswith( "#"):
        parts = [x.strip() for x in l.split(";")]
        if len( parts) != 3:
          raise ValueError( "wrong line in data file: '%s'" % l)
        name, smiles_string, ring_hash = parts
        rng = ring( name, smiles_string, ring_hash=ring_hash)
        self.rings[ ring_hash] = rng

  def _analyze_structure_dependencies( self):
    for v1 in self.structures.vertices:
      for v2 in self.structures.vertices:
        if v1 is not v2:
          sub1 = v1.value
          sub2 = v2.value
          if sub1.structure.contains_substructure( sub2.structure):
            self.structures.add_edge( v2, v1)

  def _find_head_structures( self):
    for v1 in self.structures.vertices:
      v1.properties_['in_links'] = []
      for v2 in self.structures.vertices:
        if v1 in v2.neighbors:
          v1.properties_['in_links'].append( v2.value)
    heads = []
    for v in self.structures.vertices:
      if len( v.properties_['in_links']) == 0:
        # we have got a top of the tree
        heads.append( v)
    return heads
  
  def _compute_search_trees( self):
    heads = self._find_head_structures()
    assert heads
    for head in heads:
      d = self.structures.mark_vertices_with_distance_from( head)
      dv = [(v.properties_['d'],v) for v in self.structures.vertices if 'd' in v.properties_]
      dv.sort( reverse=True)
      for d,v in dv:
        for d,parent in dv:
          if parent.value in v.properties_['in_links']:
            if v.value not in parent.value.children:
              # this ensures only one copy of a substructure in a tree
              parent.value.children.append( v.value)
            break
    return [h.value for h in heads]

  def which_substructure_is_more_specific( self, s1, s2):
    """returns 1 if s1 is more specific,
               2 if s2 is more specific,
               3 if both are the same,
               4 if the relationship cannot be established
               (there is not connection between s1 and s2 in the
               substructure interdependency graph - this signifies
               a problem, probably in the matching algorithm
               """
    if s1 is s2:
      return 3
    # just find vertices corresponding to s1 and s2
    v1 = None
    v2 = None
    for v in self.structures.vertices:
      if v.value == s1:
        v1 = v
      elif v.value == s2:
        v2 = v
    assert v1
    assert v2
    p1 = self.structures.path_exists( v1, v2)
    p2 = self.structures.path_exists( v2, v1)
    if (p1 and p2) or (not p1 and not p2):
      # both paths exist or none does - this should not happen
      return 4
    if p1:
      return 2
    if p2:
      return 1

  def find_substructures_in_mol( self, mol):
    # get the hits
    hits2 = []
    for v in self.structures.vertices:
      struct = v.value
      ms = struct.find_matches( mol)
      hits2 += ms
    # weed out the hits that match inside a ring
    ring_hits = self.find_rings_in_mol( mol)
    hits = []
    for hit in hits2:
      keep = True
      for rhit in ring_hits:
        if set(rhit.get_significant_atoms()) & set(hit.get_significant_atoms()):
          keep = False
          break
      if keep:
        hits.append( hit)
    # weed out overlapping hits - leave only the most significant ones
    hit_num = len( hits)
    to_delete = True
    while to_delete:
      to_delete = []
      for i,hit1 in enumerate( hits):
        for hit2 in hits[i+1:]:
          if hit1 is not hit2:
            if set(hit1.get_significant_atoms()) & set(hit2.get_significant_atoms()):
              winner = self.which_substructure_is_more_specific( hit1.substructure, hit2.substructure)
              if winner == 1:
                to_delete.append( hit2)
              elif winner == 2:
                to_delete.append( hit1)
              elif winner == 3:
                pass # we preserve both hits
              else:
                raise ValueError( "Relationship between competing fragments could not be established,\nthere is probaly an error in the substructure matching code.")
        if to_delete:
          break
      hits = [hit for hit in hits if not hit in to_delete]
    return ring_hits + hits

  def find_rings_in_mol( self, mol, detect_aromatic=True):
    hits = []
    if detect_aromatic:
      # this is necessary to correctly process fused aromatic rings
      # with improperly localized bonds
      mol.mark_aromatic_bonds()
    for ering in mol.get_smallest_independent_cycles_e():
      vring = mol.edge_subgraph_to_vertex_subgraph( ering)
      ring_mol = mol.get_new_induced_subgraph( vring, ering)
      # here we need to take care of aromatic bonds
      # this is needed to properly match naphthalene, etc.
      if detect_aromatic:
        for e in ring_mol.edges:
          if e.aromatic:
            e.order = 4
        ring_mol.localize_aromatic_bonds()
      ring_hash = ring_mol.get_structure_hash()
      if ring_hash in self.rings:
        hit = ring_match( vring, self.rings[ ring_hash])
      else:
        ring_obj = ring( None, smiles.mol_to_text( ring_mol), ring_hash=ring_hash)
        self.rings[ ring_hash] = ring_obj
        hit = ring_match( vring, ring_obj)
      hits.append( hit)
    return hits

  @classmethod
  def _data_files_to_python_module( self, structure_file=None, ring_file=None):
    out = file( "subsearch_data.py", "w")
    f = file( structure_file or self.substructure_def_file, "r")
    print("""#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------

""", file=out)
    print("## automatically generated file - may be overwritten at any time", file=out)
    print("structures = [", file=out)
    for line in f:
      l = line.strip()
      if l and not l.startswith( "#"):
        parts = [x.strip() for x in l.split(";")]
        if len( parts) < 3:
          print("Invalid line in src file:", line[:-1], file=sys.stderr)
        elif len( parts) == 3:
          parts.append( "")
        to_ignore = map( int, filter( None, parts[3].split(",")))
        parts[3] = to_ignore
        if not parts[1]:
          parts[1] = parts[0]
        print(tuple(parts), ",", file=out)
    print("]", file=out)
    f.close()
    f = file( ring_file or self.ring_def_file, "r")
    print("rings = [", file=out)
    for line in f:
      l = line.strip()
      if l and not l.startswith( "#"):
        parts = [x.strip() for x in l.split(";")]
        print(tuple(parts), ",", file=out)
    print("]", file=out)
    f.close()
    out.close()


class substructure( object):

  def __init__( self, name, compound_type, smiles="", atoms_to_ignore=None):
    self.name = name
    self.compound_type = compound_type
    self.structure = None
    if smiles:
      self.read_smiles( smiles, atoms_to_ignore=atoms_to_ignore)
    self.children = []

  def __str__( self):
    return "<Substructure: %s: %s>" % (self.name, self.smiles_string)

  def read_smiles( self, smiles_string, atoms_to_ignore=None):
    self.smiles_string = smiles_string.strip()
    self.structure = smiles.text_to_mol( smiles_string, calc_coords=False)
    if atoms_to_ignore:
      self.atoms_to_ignore = [self.structure.vertices[x-1] for x in atoms_to_ignore]
    else:
      self.atoms_to_ignore = []

  def find_matches( self, mol):
    ret = []
    ms = list( mol.select_matching_substructures( self.structure, implicit_freesites=True, auto_cleanup=False))
    for atoms in ms:
      num = min( atoms[0].properties_['subsearch'].keys())
      atoms_in_fragment = [a.properties_['subsearch'][num] for a in atoms]
      ret.append( substructure_match( atoms, atoms_in_fragment, self))
    mol.clean_after_search( self.structure)
    return ret


class ring( object):

  def __init__( self, name, smiles, ring_hash=None):
    self.name = name
    self.compound_type = "ring"
    self.smiles_string = smiles.strip()
    if ring_hash:
      self.ring_hash = ring_hash
    else:
      self.ring_hash = None # read smiles and generate hash here?

  def __str__( self):
    return "<Ring: %s: %s>" % (self.name, self.smiles_string)


class substructure_match( object):

  def __init__( self, atoms_found, atoms_searched, substruct):
    """atoms_found are atoms in the molecule we searched in,
    atoms_searched are atoms in the fragment we used for search,
     (atoms_searched and atoms_found are guaranteed to have the same order, thus
      allowing matching between the two structures),
    substruct is the substructure instance"""
    self.substructure = substruct
    self.atoms_found = atoms_found
    self.atoms_searched = atoms_searched

  def __str__( self):
    return "<Match of %s with %d atoms (%d significant)>" % (self.substructure, len( self.atoms_found), len( self.get_significant_atoms()))

  def get_significant_atoms( self):
    ret = []
    for i,af in enumerate( self.atoms_found):
      atms = self.atoms_searched[i]
      if atms not in self.substructure.atoms_to_ignore:
        ret.append( af)
    return ret


class ring_match( object):

  def __init__( self, atoms_found, ring_obj):
    self.substructure = ring_obj
    self.atoms_found = atoms_found

  def __str__( self):
    return "<RingMatch of %s>" % (self.substructure)

  def get_significant_atoms( self):
    return self.atoms_found


  
if __name__ == "__main__":
  # update the subsearch_data.py module
  substructure_search_manager._data_files_to_python_module()
  import time
  t = time.time()
  ssm = substructure_search_manager()
  #ssm._read_structure_file()
  #ssm._read_ring_file()
  
  print("Read_structure_file: %.1fms" % (1000*(time.time() - t)))
  t = time.time()
  
  print(ssm.structures)
  #print ssm.structures.is_connected()
  dump = ssm.structures.get_graphviz_text_dump()
  f = file( "dump.dot", "w")
  f.write( dump)
  f.close()

  print("Graphviz dump: %.1fms" % (1000*(time.time() - t)))
  t = time.time()

  def print_tree( x, l):
    print(l*" ", x)  #, x.children
    for ch in x.children:
      print_tree( ch, l+2)
  
  #for tree in ssm._compute_search_trees():
  #  print_tree( tree, 0)

  #"c1cccc2c1cccc2
  #text = "COc5ccc4c2sc(cc2nc4c5)-c(cc1nc3c6)sc1c3ccc6OC"
  text = 'C(=O)OCC'
  #text2 = 'C(=O)[O-]'
  mol = smiles.text_to_mol( text, calc_coords=False)
  #m2 = smiles.text_to_mol( text2, calc_coords=False)
  #print "XXXX", mol.contains_substructure( m2)

  subs = ssm.find_substructures_in_mol( mol)
  for sub in subs:
    print(sub)

  print("Substructure search: %.1fms" % (1000*(time.time() - t)))
  t = time.time()
