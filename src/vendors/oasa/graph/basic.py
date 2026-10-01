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
#
# Last edited: $Date: 2003/12/23 18:06:15 $
#
#--------------------------------------------------------------------------


class attribute_flexible_class:
  """this class provides mechanism for addition and removing of attributes on fly"""

  def __init__( self):
    pass

  def add_attribute( self, name, value=None):
    if name not in self.__dict__:
      self.__dict__[ name] = value
      return 1
    return 0

  def del_attribute( self, name):
    if name not in self.__dict__:
      del self.__dict__[ name]
      return 1
    return 0
    
