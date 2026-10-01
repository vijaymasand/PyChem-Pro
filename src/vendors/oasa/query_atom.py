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

import sys
sys.path.append( '../')

from . import graph
from . import periodic_table as PT
from .common import is_uniquely_sorted

import copy
from warnings import warn

import re
from .atom import atom

from .chem_vertex import chem_vertex
from .oasa_exceptions import oasa_invalid_atom_symbol



class query_atom( chem_vertex):

  ## ("value","charge","x","y","z","multiplicity","valency","charge","free_sites")
  attrs_to_copy = chem_vertex.attrs_to_copy + ("symbols",)

  def __init__( self, coords=None):
    chem_vertex.__init__( self, coords=coords)
    self.symbols = set()
    self.free_sites = 0


  def matches( self, other):
    if not isinstance( other, atom):
      return False
    #if isinstance( other, query_atom):
    #  raise ValueError, "query atoms on the recieving side are not supported"
    # halogens
    for sym in self.symbols:
      if sym == other.symbol:
        return True
      elif sym == "X":
        if other.symbol in ("F","Cl","Br","I"):
          return True
      # Q (any except H,C)
      elif sym == "Q":
        if other.symbol not in "HC":
          return True
      # A (any except H)
      elif sym == "A":
        if other.symbol != "H":
          return True
      # R - anything
      elif sym == "R":
        return True

    return False


  ## PROPERTIES


  # symbol
  def _set_symbol( self, symbol):
    if symbol in PT.periodic_table.keys():
      if not "query" in PT.periodic_table[ symbol]:
        warn( "Setting normal atom symbol to a query_atom instance, do you mean it?")
      self.symbols = set( [symbol])
    else:
      self.symbols = self.parse_query_definition( symbol)

    self._symbol = symbol


  def _get_symbol( self):
    return self._symbol

  symbol = property( _get_symbol, _set_symbol, None, "atom symbol")





  # free_sites
  def _set_free_sites( self, free_sites):
    self._free_sites = free_sites

  def _get_free_sites( self):
    return self._free_sites

  free_sites = property( _get_free_sites, _set_free_sites, None, "atoms free_sites")



      
  def __str__( self):
    return "query atom '%s'" % str( self.symbol)




  def is_query_definition( self, text):
    matcher = re.compile( r"\[([A-Z][a-z]?,)*[A-Z][a-z]?\]")
    return matcher.match( text) and True or False

  is_query_definition = classmethod( is_query_definition)


  def parse_query_definition( self, text):
    if self.is_query_definition( text):
      syms = set( map( str, text[1:-1].split(",")))
      for sym in syms:
        if sym not in PT.periodic_table.keys():
          raise oasa_invalid_atom_symbol( "invalid symbol in query definition", sym)
      return syms
    else:
      raise oasa_invalid_atom_symbol( "not valid query definition", text)

  parse_query_definition = classmethod( parse_query_definition)



#print query_atom.parse_query_definition( "[Ar,Cl,X]")
