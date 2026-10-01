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

from .config import Config

import re
from .periodic_table import periodic_table
from . import coords_generator
from . import misc
from .known_groups import name_to_smiles
from . import smiles
from .oasa_exceptions import oasa_invalid_atom_symbol


class linear_formula( object):

  def __init__( self, text="", start_valency=0, end_valency=0, mol=None):
    """valency specifies the already occupied valency of the submited formula,
    it is usually used when parsing groups"""
    self.molecule = None
    if text:
      self.parse_text( text, start_valency=start_valency, end_valency=end_valency, mol=mol)


  def parse_text( self, text, start_valency=0, end_valency=0, mol=None):
    text = self.expand_abbrevs( text)
    mol = self.parse_form( text, start_valency=start_valency, mol=mol, reverse=False)
    if mol:
      # are there any atoms?
      if not mol.vertices:
        return None
      
      # now we check if the structure is complete
      for v in mol.vertices:
        if v.free_valency:
          # we permit free valency at the last atom if end_valency was given
          if not (v == self.get_last_free_atom(mol) and v.free_valency == end_valency):
            return None

      if start_valency:
        mol.remove_vertex( mol.vertices[0]) # remove the dummy

      # are there any atoms, again?
      if not mol.vertices:
        return None

      # is the molecule connected?
      if not mol.is_connected():
        return None

      # remove hydrogens, but only if the molecule is not completely made of them :)
      mol.remove_unimportant_hydrogens()

      # set some objects attributes for later reading
      self.molecule = mol
      if start_valency:
        self.first_atom = mol.vertices[0]
      else:
        self.first_atom = None
      if end_valency:
        self.last_atom = self.get_last_free_atom(mol)
      else:
        self.last_atom = None
      return mol


  def parse_form( self, text, start_valency=0, mol=None, reverse=False):
    form = text

    # the code itself
    if not mol:
      mol = Config.create_molecule()

    # create the dummy atom
    if start_valency:
      last_atom = mol.create_vertex()
      last_atom.valency = start_valency
      mol.add_vertex( last_atom)
    else:
      last_atom = None

    # check if there are branches in the formula
    if "(" not in form:
      # there are no subbranches
      chunks = re.split( "([A-Z][a-z]?[0-9]?[+-]?)", form)
      if reverse:
        chunks = reverse_chunks( chunks)
      for chunk in chunks:
        if chunk:
          atms = self.chunk_to_atoms( chunk, mol)
          if atms == None:
            return None
          last_atom = self.get_last_free_atom( mol)
          for a in atms:
            mol.add_vertex( a)
            if last_atom:
              max_val = min( last_atom.free_valency, a.free_valency, 3)
              if max_val <= 0 and last_atom.free_valency <= 0:
                if last_atom.raise_valency():
                  max_val = min( last_atom.free_valency, a.free_valency, 3)
              b = mol.create_edge()
              b.order = max_val
              mol.add_edge( last_atom, a, b)
            else:
              last_atom = a
    else:
      for chunk, count in gen_formula_fragments( form, reverse=reverse):
        if chunk:
          last_atom = self.get_last_free_atom( mol)
          do_linear = False # should we string the fragments rather than adding them all to the last atom
          if last_atom and last_atom.free_valency < count:
            do_linear = True

          for j in range( count):
            if chunk[0] == "!":
              # the form should be a smiles
              m = smiles.text_to_mol( chunk[1:], calc_coords=0)
              m.add_missing_hydrogens()
              hs = [v for v in m.vertices[0].neighbors if v.symbol == 'H']
              m.disconnect( hs[0], m.vertices[0])
              m.remove_vertex( hs[0])
              smile = True
            else:
              val = last_atom and 1 or 0
              m = self.parse_form( chunk, start_valency=val, mol=mol.create_graph())
              smile = False
            if not m:
              return None

            if not last_atom:
              # !!! this should not happen in here
              mol.insert_a_graph( m)
              # if there are multiple chunks without a previous atom, we should create a linear fragment
              #last_atom = [v for v in mol.vertices if v.free_valency > 0][-1]
            else:
              if not smile:
                m.remove_vertex( m.vertices[0]) # remove the dummy
              mol.insert_a_graph( m)
              b = mol.create_edge()
              mol.add_edge( last_atom, m.vertices[0], b)
              if do_linear:
                last_atom = m.vertices[0]

    return mol
        


  def chunk_to_atoms( self, chunk, mol):
    m = re.match( "([A-Z][a-z]?)([0-9])?([+-])?", chunk)    
    if m:
      name = m.group( 1)
      number = m.group( 2) and int( m.group(2)) or 1
      sign = m.group( 3) and (int( m.group(3)+'1')) or 0
      ret = []
      for i in range( number):
        v = mol.create_vertex()
        try:
          v.symbol = name
        except oasa_invalid_atom_symbol:
          return None
        v.charge = sign
        ret.append( v)
      return ret


  def get_last_free_atom( self, mol):
    # check if there is something with a free valency
    if not mol.vertices:
      return None
    atoms = [o for o in misc.reverse( mol.vertices)]
    for a in atoms:
      if a.free_valency > 0:
        return a
    # if its not the case
    for i, a in enumerate( atoms[0:-1]):
      b = a.get_edge_leading_to( atoms[i+1])
      if b and b.order > 1:
        b.order -= 1
        return a
    # well, we cannot do anything else
    return atoms[-1]


  def expand_abbrevs( self, text):
    # at first sort the text according to length, so that the longest are expanded first
    # (MMTr and not Tr)
    keys = [(len( k), k) for k in name_to_smiles.keys()]
    keys.sort()
    keys.reverse()
    for l, key in keys:
      val = name_to_smiles[ key]
      text = text.replace( key, "(!%s)" % val)
    return text


def reverse_chunks( chunks):
  chunks.reverse()
  for i in range( 0, len( chunks), 2):
    chunks[i], chunks[i+1] = chunks[i+1],chunks[i]
  return chunks


def gen_formula_fragments( formula, reverse=False):
  chunks = list( gen_formula_fragments_helper( formula))
  if reverse:
    chunks.reverse()
  i = 0
  while i < len( chunks):
    chunk, brack = chunks[ i]
    if brack and i < len( chunks) - 1:
      next, nbrack = chunks[i+1]
      if not nbrack:
        count, rest = split_number_and_text( next)
        if count == None:
          yield chunk, 1
        else:
          chunks[i+1] = (rest, nbrack)
          if not rest:
            i += 1
          yield chunk, count
      else:
        yield chunk, 1
    else:
      yield chunk, 1
    i += 1


def gen_formula_fragments_helper( formula):
  opened_brackets = 0
  to_ret = []
  for ch in formula:
    if ch not in "()":
      to_ret.append( ch)
    elif ch == "(":
      if opened_brackets == 0:
        if to_ret:
          yield ''.join( to_ret), False
          to_ret = []
      else:
        to_ret.append( "(")
      opened_brackets += 1
    elif ch == ")":
      opened_brackets -= 1
      if opened_brackets < 0:
        opened_brackets = 0
      if opened_brackets == 0:
        if to_ret:
          yield ''.join( to_ret), True
          to_ret = []
      else:
        to_ret.append( ")")
  if to_ret:
    yield ''.join( to_ret), False



def split_number_and_text( txt):
  last = None
  for i in range( len( txt)):
    try:
      last = int( txt[0:i+1])
    except ValueError:
      return last, txt[i:]
  return last, ""
    


def reverse_formula( text):
  all_chunks = []
  form = text
  if "(" not in form:
    chunks = re.split( "([A-Z][a-z]?[0-9]?[+-]?)", form)
    for chunk in chunks:
      if chunk:
        all_chunks.append( chunk)
  else:
    for chunk, count in gen_formula_fragments_helper( form):
      if count:
        all_chunks.append(")")
      if chunk:
        all_chunks.extend( reverse_formula( chunk))
      if count:
        all_chunks.append("(")
#
#        if count:
#          all_chunks.append( str(count))

  print(all_chunks)
  #all_chunks.reverse()
  return all_chunks
  
  
if __name__ == "__main__":
  form = 'CH3(CH2)7'
  #form = 'CH3((CH2)2)2O'
  form = "1) Bu"

  a = linear_formula( form, start_valency=0, end_valency=1)
  if a.molecule:
    m = a.molecule
    coords_generator.calculate_coords( m)
    coords_generator.show_mol( m)
  else:
    print("nothing")

