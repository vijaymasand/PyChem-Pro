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


from math import pi, sqrt, sin, cos, atan2
from . import geometry
from . import misc

import warnings
#import graph


def sum_of_ring_internal_angles( size):
  return (size-2)*180

def optimal_ring_iternal_angle( size, angle_shift=0):
  return sum_of_ring_internal_angles( size)/size


def gen_variations(items, n):
    if n==0:
      yield []
    else:
      for i in range( len(items)-n+1):
        for v in gen_variations(items[i+1:],n-1):
          yield [items[i]]+v


def deg_to_rad( deg):
  return pi*deg/180

def rad_to_deg( rad):
  return 180*rad/pi




class coords_optimizer:

  def __init__( self):
    self.cutoff_rms_grad = 1e-5
    self.cutoff_max_grad = 3e-5
    self.max_iter_number = 1000


  def optimize_coords( self, mol, bond_length=1, callback=None):
    """callback may be used to obtain information about the running optimization,
    it is called after each step with three parameters - step number, RMS grad and maxgrad"""
    self.mol = mol

    if bond_length < 0:
      bls = []
      for b in mol.edges:
        a1, a2 = b.vertices
        bls.append( sqrt( (a1.x-a2.x)**2 + (a1.y-a2.y)**2))
      self.bond_length = sum( bls) / len( bls)
    else:
      self.bond_length = bond_length

    sumdd, max_grad = self.step()
    rms_grad = sumdd / len( self.mol.vertices)
    i = 0
    ok = True
    self.end_reason = "converged"
    
    while rms_grad > (self.cutoff_rms_grad / self.bond_length) or max_grad > (self.cutoff_max_grad / self.bond_length):
      i += 1
      if callback:
        callback( i, rms_grad, max_grad)

      sumdd, max_grad = self.step()
      rms_grad = sumdd / len( self.mol.vertices)

      if rms_grad > 1000:
        self.end_reason = "gradien too large"
        ok = False
        break
      if i > self.max_iter_number:
        self.end_reason = "max_iter_number iterations without converged result"
        ok = False
        break

    self.i = i
    self.rms_grad = rms_grad
    self.max_grad = max_grad
    return ok


  def step( self, callback=None):

    for v in self.mol.vertices:
      v.properties_['dx'] = 0
      v.properties_['dy'] = 0

    for v in self.mol.vertices:
      self.set_gradient( v)

    dd = 0
    mdd = 0
    for v in self.mol.vertices:
      dx = v.properties_['dx']
      dy = v.properties_['dy']
      v.x += dx
      v.y += dy
      d = dx**2 + dy**2
      dd += d
      if d > mdd:
        mdd = d

    return sqrt( dd), sqrt( mdd)



  def set_gradient( self, vertex):
    dx, dy = 0, 0
    for v in vertex.neighbors:
      gx, gy = self.get_length_gradient( vertex, v)
      dx -= gx/2.0
      dy -= gy/2.0
    vertex.properties_['dx'] += dx
    vertex.properties_['dy'] += dy
    
    dx, dy = 0, 0
    vs = list( vertex.neighbors)
    variations = list( gen_variations( vs, 2))
    if variations:
      for var in variations:
        v1, v2 = var
        gx1, gy1, gx2, gy2 = self.get_angle_gradient( deg_to_rad( 120), vertex, v1, v2)
        v1.properties_['dx'] -= gx1/4.0
        v1.properties_['dy'] -= gy1/4.0
        v2.properties_['dx'] -= gx2/4.0
        v2.properties_['dy'] -= gy2/4.0


  def get_length_gradient( self, v1, v2, opt_length=0):
    if not opt_length:
      opt_length = self.bond_length
    dx = v1.x - v2.x
    dy = v1.y - v2.y
    dist = sqrt( dx**2 + dy**2)
    dd = dist - opt_length
    gx = (dx / dist) * dd
    gy = (dy / dist) * dd
    return gx, gy


  def get_angle_gradient( self, opt_angle, refv, v1, v2):
    expl = sqrt( 3) * self.bond_length
    gx1, gy1 = self.get_length_gradient( v1, v2, opt_length=expl)
    gx2 = -gx1
    gy2 = -gy1
    return gx1, gy1, gx2, gy2


  def get_angle_gradient2( self, opt_angle, refv, v1, v2):
    ang1 = geometry.clockwise_angle_from_east( v1.x-refv.x, v1.y-refv.y)
    ang2 = geometry.clockwise_angle_from_east( v2.x-refv.x, v2.y-refv.y)
    ang = ang1 - ang2
    sign = -geometry.on_which_side_is_point( (refv.x, refv.y, v1.x, v1.y), (v2.x, v2.y))
    
    while ang < 0:
      ang += 2*pi
    if ang > pi:
      ang = 2*pi - ang
      sign *= -1

    dang = (ang - opt_angle) / 2

    gx1 = (cos( ang1 - sign*dang) - cos( ang1)) * sqrt( (refv.x-v1.x)**2 + (refv.y-v1.y)**2) #self.bond_length #(v1.x-refv.x)
    gy1 = (sin( ang1 - sign*dang) - sin( ang1)) * sqrt( (refv.x-v1.x)**2 + (refv.y-v1.y)**2) #self.bond_length #(v1.y-refv.y)

    gx2 = (cos( ang2 + sign*dang) - cos( ang2)) * sqrt( (refv.x-v2.x)**2 + (refv.y-v2.y)**2) #self.bond_length #(v2.x-refv.x)
    gy2 = (sin( ang2 + sign*dang) - sin( ang2)) * sqrt( (refv.x-v2.x)**2 + (refv.y-v2.y)**2) #self.bond_length #(v2.y-refv.y)


    # control
    ang1 = geometry.clockwise_angle_from_east( v1.x-refv.x+gx1, v1.y-refv.y+gy1)
    ang2 = geometry.clockwise_angle_from_east( v2.x-refv.x+gx2, v2.y-refv.y+gy2)
    ang = ang1 - ang2
    sign = geometry.on_which_side_is_point( (refv.x, refv.y, v1.x, v1.y), (v2.x, v2.y))
    
    while ang < 0:
      ang += 2*pi
      #sign *= -1
    if ang > pi:
      ang = 2*pi - ang
      #sign *= -1

    dang2 = (ang - opt_angle) / 2

    if abs( dang2) > abs( dang):
      print("fuck", rad_to_deg( ang1), rad_to_deg( ang2), rad_to_deg( dang))
    else:
      #rint "good", rad_to_deg( ang1), rad_to_deg( ang2), rad_to_deg( dang)
      pass

    return gx1, gy1, gx2, gy2




if __name__ == "__main__":

  from . import smiles
  import time

  mol = smiles.text_to_mol( 'CC(C)C(C)CCCC')

  print([(v.x,v.y) for v in mol.vertices])

  opt = coords_optimizer()
  t = time.time()
  ok = opt.optimize_coords( mol, bond_length=0.9)
  print("converged:", ok)
  print("%d iterations, RMS grad %f, max grad %f" % (opt.i, opt.rms_grad, opt.max_grad))
  print("time %.1f ms" % ((time.time() - t) * 1000))


  print([(v.x,v.y) for v in mol.vertices])

  from . import coords_generator

  coords_generator.show_mol( mol)
