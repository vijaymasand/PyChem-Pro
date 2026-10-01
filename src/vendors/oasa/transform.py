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
#
#
#
#--------------------------------------------------------------------------

"""this module provides basic coordinate transformations based on matrix algebra"""


from math import cos, sin, pi
from . import geometry

class transform:
  """this class provides basic higher-level interface for coordinate transforms"""
  def __init__( self, mat = None):
    if not mat:
      self.mat = matrix( [[1,0,0],[0,1,0],[0,0,1]])
    else:
      self.mat = matrix( mat)

  def transform_xy( self, x, y):
    x1, y1, one = self.mat.get_multiplied2( [[x],[y],[1]])
    return x1[0], y1[0]

  def transform_xyz( self, x, y, z):
    """used for compatibility with transform3d, it does not work with
    z-coord at all"""
    x,y = self.transform_xy( x, y)
    return x,y,z

  def transform_4( self, list):
    "for items that have 4 coordinates (lines, rectangles) this transforms 'em all"
    x1, y1, one = self.mat.get_multiplied2( [[list[0]],[list[1]],[1]])
    x2, y2, one = self.mat.get_multiplied2( [[list[2]],[list[3]],[1]])
    return (x1[0], y1[0], x2[0], y2[0])

  def transform_list( self, l):
    ret = []
    for line in l:
      ret.append( self.transform_xy( line[0], line[1]))
    return ret

  def transform_xy_flat_list( self, coords):
    """transforms a list that cointains alternating x,y values (not list of pairs
    as self.transform_list)"""
    ret = []
    for j in range( 0, len( coords), 2):
      x, y = self.transform_xy( coords[j], coords[j+1])
      ret.append( x)
      ret.append( y)
    return ret


  def map_transform( self, l):
    "this method can be used with map() function"
    return self.transform_xy( l[0], l[1])

  def set_move( self, dx, dy):
    "add an moving step to transformation matrix"
    self.mat = matrix( mat=self.mat.get_multiplied( [[1,0,dx],[0,1,dy],[0,0,1]]))

  def set_rotation( self, angle):
    "add an rotation step to transformation matrix"
    self.mat = matrix( mat=self.mat.get_multiplied( [[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]]))

  def set_scaling_xy( self, mx, my):
    "add an scaling step to transformation matrix"
    self.mat = matrix( mat=self.mat.get_multiplied( [[mx,0,0],[0,my,0],[0,0,1]]))

  def set_scaling( self, scale):
    "add an scaling step to transformation matrix, same scaling for both dimensions"
    self.mat = matrix( mat=self.mat.get_multiplied( [[scale,0,0],[0,scale,0],[0,0,1]]))

  def get_scaling( self):
    """computes the value of scaling from self"""
    x01, y01, x02, y02 = [0, 0, 100, 100]
    x11, y11, x12, y12 = self.transform_4( (x01, y01, x02, y02))
    l1 = geometry.point_distance( x01, y01, x02, y02)
    l2 = geometry.point_distance( x11, y11, x12, y12)
    return l2/l1


  def get_scaling_xy( self):
    """computes the scaling for 'x' and 'y' separately"""
    x01, y01, x02, y02 = [0, 0, 100, 100]
    x11, y11, x12, y12 = self.transform_4( (x01, y01, x02, y02))
    return (x12-x11)/(x02-x01), (y12-y11)/(y02-y01)
 

class matrix:
  "this class provides common operations for matrix of 3x3 elements"

  def __init__( self, mat):
    self.mat = mat

  def get_transposed( self):
    m = self.mat
    ret = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range( len( m)):
      for j in range( len( m)):
        ret[j][i] = m[i][j]
    return ret

  def get_inverse( self):
    pass

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
          m[i].append( self.mat[i][0]*multi[0][j] + self.mat[i][1]*multi[1][j] + self.mat[i][2]*multi[2][j])
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
          m[i].append( multi[i][0]*self.mat[0][j] + multi[i][1]*self.mat[1][j] + multi[i][2]*self.mat[2][j])
    return  m

  def get_determinant( self):
    m = self.mat
    return (((m[0][0] * m[1][1] * m[2][2]) + (m[0][1] * m[1][2] * m[2][0]) + (m[0][2] * m[1][0] * m[2][1])) -
            ((m[2][1] * m[1][2] * m[0][0]) + (m[2][2] * m[1][0] * m[0][1]) + (m[2][0] * m[1][1] * m[0][2])))
                                                                                         
