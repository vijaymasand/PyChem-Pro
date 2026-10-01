#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2009 Beda Kosata <beda@zirael.org>

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


"""this module provides basic coordinate transformations in 3D based on matrix algebra"""

from math import cos, sin, pi

class transform3d:
  """this class provides basic higher-level interface for coordinate transforms"""
  def __init__( self, mat = None):
    if not mat:
      self.mat = matrix( [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    else:
      self.mat = matrix( mat)

  def transform_xyz( self, x, y, z):
    x1, y1, z1, one = self.mat.get_multiplied2( [[x],[y],[z],[1]])
    return x1[0], y1[0], z1[0]

  def transform_xy( self, x, y):
    """is used only for compatibility with 2D transform, assumes z=0,
    should be used only in cases where it is sure that starting and
    result positions have z0 (pi rotation around one axis etc.)"""
    x0, y0, z0 = self.transform_xyz( x, y, 0)
    return x0, y0

  def transform_xy_flat_list( self, coords):
    """ONLY X,Y TRANSFORM, for compatibility with 2D transform!!!,
    transforms a list that cointains alternating x,y values (not list of pairs
    as self.transform_list)"""
    ret = []
    for j in range( 0, len( coords), 2):
      x, y = self.transform_xy( coords[j], coords[j+1])
      ret.append( x)
      ret.append( y)
    return ret

  def transform_xyz_flat_list( self, coords):
    ret = []
    for j in range( 0, len( coords), 3):
      ret += self.transform_xyz( coords[j], coords[j+1], coords[j+2])
    return ret

  def transform_list( self, l):
    ret = []
    for line in l:
      ret.append( self.transform_xyz( line[0], line[1], line[2]))
    return ret

  def map_transform( self, l):
    "this method can be used with map() function"
    return self.transform_xyz( l[0], l[1], l[2])

  def set_move( self, dx, dy, dz):
    "add an moving step to transformation matrix"
    self.mat = matrix( mat=self.mat.get_multiplied( [[1,0,0,dx],[0,1,0,dy],[0,0,1,dz],[0,0,0,1]]))

  def set_rotation( self, xa, ya, za):
    "add an rotation step to transformation matrix"
    self.set_rotation_x( xa)
    self.set_rotation_y( ya)
    self.set_rotation_z( za)
    
  def set_rotation_x( self, xa):
    mat = [[1,0,0,0],
           [0, cos(xa), sin(xa), 0],
           [0, -sin(xa), cos(xa), 0],
           [0,0,0,1]]
    self.mat = matrix( mat=self.mat.get_multiplied( mat))

  def set_rotation_y( self, ya):
    mat = [[cos(ya), 0, -sin(ya), 0],
           [0, 1, 0, 0],
           [sin(ya), 0, cos(ya), 0],
           [0,0,0,1]]
    self.mat = matrix( mat=self.mat.get_multiplied( mat))

  def set_rotation_z( self, za):
    mat = [[cos(za), sin(za), 0, 0],
           [-sin(za), cos(za), 0, 0],
           [0,0,1,0],
           [0,0,0,1]]
    self.mat = matrix( mat=self.mat.get_multiplied( mat))

  def set_scaling_xyz( self, mx, my, mz):
    "add an scaling step to transformation matrix"
    self.mat = matrix( mat=self.mat.get_multiplied( [[mx,0,0,0],[0,my,0,0],[0,0,mz,0],[0,0,0,1]]))

  def set_scaling( self, scale):
    "add an scaling step to transformation matrix, same scaling for both dimensions"
    self.set_scaling_xyz( scale, scale, scale)

  def get_inverse( self):
    return transform3d( mat=self.mat.get_inverse())


class matrix:
  "this class provides common operations for matrix of 3x3 elements"

  def __init__( self, mat):
    self.mat = mat

  def get_transposed( self):
    m = self.mat
    ret = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    for i in range( len( m)):
      for j in range( len( m)):
        ret[j][i] = m[i][j]
    return ret

  def get_inverse( self):
    def _part( a, b):
      _ret = [[0,0,0],[0,0,0],[0,0,0]]
      for i in range( len( m)):
        if i == a:
          continue
        elif i > a:
          i2 = i - 1
        else:
          i2 = i
        for j in range( len( m)):
          if j == b:
            continue
          elif j > b:
            j2 = j - 1
          else:
            j2 = j
          _ret[i2][j2] = m[i][j]
      return _ret
      
    m = self.mat
    inv = matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    det = self.get_determinant()
    for i in range( len( m)):
      for j in range( len( m)):
        part = _part( i, j)
        part_det = matrix._determinant_3( part)
        sign = (i+j)%2 and -1.0 or 1.0
        inv.mat[i][j] = sign * part_det / det
    return inv.get_transposed()


  def get_multiplied2( self, multi):
    "returns matrix multiplied by multi"
    m = []
    if type( multi) != type([]):
      for i in range( len( self.mat)):
        m.append([])
        for j in range( len( self.mat)):
          m[i].append( multi*self.mat[i][j])
    else:
      for i in range( len( self.mat)):
        m.append([])
        for j in range( len( multi[0])):
          m[i].append( self.mat[i][0]*multi[0][j] + self.mat[i][1]*multi[1][j] + self.mat[i][2]*multi[2][j] + self.mat[i][3]*multi[3][j])
    return  m

  def get_multiplied( self, multi):
    "returns matrix multiplied by multi"
    m = []
    if type( multi) != type([]):
      for i in range( len( self.mat)):
        m.append([])
        for j in range( len( self.mat)):
          m[i].append( multi*self.mat[i][j])
    else:
      for i in range( len( multi)):
        m.append([])
        for j in range( len( self.mat)):
          m[i].append( multi[i][0]*self.mat[0][j] + multi[i][1]*self.mat[1][j] + multi[i][2]*self.mat[2][j] +  multi[i][3]*self.mat[3][j])
    return  m

  def get_determinant( self):
    _d3 = matrix._determinant_3
    m = self.mat
    a = m[0][0] * _d3([[m[1][1],m[1][2],m[1][3]],[m[2][1],m[2][2],m[2][3]],[m[3][1],m[3][2],m[3][3]]])
    b = m[0][1] * _d3([[m[1][0],m[1][2],m[1][3]],[m[2][0],m[2][2],m[2][3]],[m[3][0],m[3][2],m[3][3]]])
    c = m[0][2] * _d3([[m[1][0],m[1][1],m[1][3]],[m[2][0],m[2][1],m[2][3]],[m[3][0],m[3][1],m[3][3]]])
    d = m[0][3] * _d3([[m[1][0],m[1][1],m[1][2]],[m[2][0],m[2][1],m[2][2]],[m[3][0],m[3][1],m[3][2]]])
    return a-b+c-d


  @staticmethod
  def _determinant_3( _m):
      return (((_m[0][0] * _m[1][1] * _m[2][2]) + (_m[0][1] * _m[1][2] * _m[2][0]) + (_m[0][2] * _m[1][0] * _m[2][1])) - ((_m[2][1] * _m[1][2] * _m[0][0]) + (_m[2][2] * _m[1][0] * _m[0][1]) + (_m[2][0] * _m[1][1] * _m[0][2])))

