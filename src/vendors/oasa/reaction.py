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

from .molecule import molecule


class reaction( object):
  """Reaction representation"""


  def __init__( self, reactants=None, products=None, reagents=None):
    self.reactants = reactants or []
    self.products = products or []
    self.reagents = reagents or []



class reaction_component( object):
  """Represents one component of a reaction"""

  def __init__( self, mol=None, stoichiometry=1):
    self.stoichiometry = stoichiometry
    self.molecule = mol

  def _set_molecule( self, mol):
    assert isinstance( mol, molecule)
    self._molecule = mol

  def _get_molecule( self):
    return self._molecule

  molecule = property( _get_molecule, _set_molecule)

  def _set_stoichiometry( self, stoich):
    assert type( stoich) in (int, float)
    self._stoichiometry = stoich

  def _get_stoichiometry( self):
    return self._stoichiometry

  stoichiometry = property( _get_stoichiometry, _set_stoichiometry)

  def __str__( self):
    return "%s * (%s)" % (self.stoichiometry, self.molecule)
