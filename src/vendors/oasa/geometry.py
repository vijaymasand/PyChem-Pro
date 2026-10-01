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


"""support module for some geometric mesurements ( geometric tramforms are in transform.py)"""


from math import sqrt, atan2, pi, cos, sin
from .misc import signum, normalize_coords
import operator



def find_parallel( x1, y1, x2, y2, d):
  "returns tuple of coordinates for parallel abscissa in distance d"
  # following is here to ensure that signum of "d" clearly determines
  # the side of line on whitch the parallel is drawn  
  if round( y2, 3) -round( y1, 3) != 0:
    if y2 < y1:
      d = -d 
    k = -(x2-x1)/(y2-y1)
    x = ( d + sqrt( k**2 +1)*x1)/ sqrt( k**2 +1)
    y = y1 + k*( x -x1)
    x0 = ( d + sqrt( k**2 +1)*x2)/ sqrt( k**2 +1)
    y0 = y2 + k*( x0 -x2)
  else:
    if x1 > x2:
      d = -d
    x, x0 = x1, x2
    y = y1 - d
    y0 = y2 - d
  return (x, y, x0, y0)

def get_parallel_distance( l1, l2):
  "returns distance of two parallels - does not check whether the two are really parallels"
  x1a, y1a, x2a, y2a = l1
  x1b, y1b, x2b, y2b = l2
  if x1a == x2a:
    return y1a-y2a
  elif y1a == y2a:
    return x1a-x2a
  else:
    ka = (y2a - y1a)/(x2a - x1a)
    qa = y1a - ka*x1a
    kb = ka
    qb = y1b - kb*x1b
    k_rev = -1/ka
    q_rev = y1a - k_rev*x1a
    x0 = (qb - q_rev)/(k_rev - kb)
    y0 = k_rev*x0 + q_rev
    d = sqrt( (x1a-x0)**2 + (y1a-y0)**2)
    if qa < qb:
      d = -d
    return d
  

def get_parallel_signum( l1, l2):
  x1a, y1a, x2a, y2a = l1
  x1b, y1b, x2b, y2b = l2
  if x1a == x2a:
    return signum( -x1a + x1b)
  elif y1a == y2a:
    return signum( -y1a+y1b)
  else:
    return signum( -x2a+x2b)

def on_which_side_is_point( line, point, threshold=0):
  """tells whether a point is on one side of a line or on the other (1,0,-1) - 0 is for point on line.
  line is given as sequence of four coordinates, point as sequence of two coords,
  threshold means what smallest angle is considered to still be on the line"""
  x1, y1, x2, y2 = line
  x, y = point
  a = atan2( y-y1, x-x1)
  b = atan2( y2-y1, x2-x1)
  if a*b < 0 and abs(a-b) > pi:
    if a < 0:
      a += 2*pi
    else:
      b += 2*pi
  if abs( a-b) <= threshold or abs( abs( a-b) -pi) <= threshold:
    return 0    
  elif a-b < 0:
    return 1
  else:
    return -1

def point_on_circle( center_x, center_y, radius, direction = (), resolution = 15):
  """finds point on circle in direction of (dx, dy), optionaly rounds the angle
  according to resolution"""
  dx, dy = direction
  if resolution:
    angle = round( atan2( dy, dx)/(pi*resolution/180.0))*(pi*resolution/180.0)
  else:
    angle = atan2( dy, dx)
  x = center_x + round( cos( angle) *radius, 2)
  y = center_y + round( sin( angle) *radius, 2)
  return x,y


def angle_between_lines( l1, l2):
  """returns angle between two lines"""
  a = clockwise_angle_from_east( l1[2]-l1[0], l1[3]-l1[1]) - clockwise_angle_from_east( l2[2]-l2[0], l2[3]-l2[1])
  if a > pi:
    return a - 2*pi
  elif a < -pi:
    return a + 2*pi
  else:
    return a


def clockwise_angle_from_east( dx, dy):
  """returns the angle in clockwise direction between the center-east line and direction"""
  angle = atan2( dy, dx)
  if angle < 0:
    angle = 2*pi + angle
  return angle



def intersection_of_line_and_rect( line, rect, round_edges=0):
  """finds a point where a line and a rectangle intersect,
  both are given as lists of len == 4"""
  lx0, ly0, lx1, ly1 = map( float, line)
  rx0, ry0, rx1, ry1 = map( float, normalize_coords( rect))

  # find which end of line is in the rect and reverse the line if needed
  if (lx0 > rx0) and (lx0 < rx1) and (ly0 > ry0) and (ly0 < ry1):
    lx0, lx1 = lx1, lx0
    ly0, ly1 = ly1, ly0

  # the computation itself
  ldx = lx1 - lx0
  ldy = ly1 - ly0

  if abs( ldx) > 0.0001:
    # we calculate using y = f(x)
    k = ldy/ldx
    q = ly0 - k*lx0
    if ldx < 0:
      xx = rx1
    else:
      xx = rx0
    xy = k*xx + q
    # the result must be in the rectangle boundaries
    # but sometimes is not because rounding problems
    if not ry0 <= xy <= ry1:
      xx = lx0
      xy = ly0
  else:
    xx = lx0
    xy = ly0
    
  if abs( ldy) > 0.0001:
    # we calculate using x = f(y)
    k = ldx/ldy
    q = lx0 - k*ly0
    if ldy < 0:
      yy = ry1
    else:
      yy = ry0
    yx = k*yy + q
    # the result must be in the rectangle boundaries
    # but sometimes is not because rounding problems
    if not rx0 <= yx <= rx1:
      yy = ly0
      yx = lx0
  else:
    yy = ly0
    yx = lx0
    
  if point_distance( lx0, ly0, xx, xy) < point_distance( lx0, ly0, yx, yy):
    return (yx, yy)
  else:
    return (xx, xy)



def point_distance( x1, y1, x2, y2):
  return sqrt( (x2-x1)**2 + (y2-y1)**2)


def line_length( _args_):
  x1,y1,x2,y2 = _args_
  return point_distance( x1, y1, x2, y2)


def rectangle_intersection( rect1, rect2):
  pass


def do_rectangles_intersect( rect1, rect2):
  xs1 = [rect1[0], rect1[2]]
  xs2 = [rect2[0], rect2[2]]  
  xs = xs1 + xs2
  ys1 = [rect1[1], rect1[3]]
  ys2 = [rect2[1], rect2[3]]
  ys = ys1 + ys2
  
  dx1 = abs( xs1[0] - xs1[1])
  dx2 = abs( xs2[0] - xs2[1])
  dx = max( xs) - min( xs)

  dy1 = abs( ys1[0] - ys1[1])
  dy2 = abs( ys2[0] - ys2[1])
  dy = max( ys) - min( ys)

  if dx1+dx2 > dx and dy1+dy2 > dy:
    return True
  else:
    return False


def is_point_inside_polygon( point, polygon):
  """point is a tuple or list of length == 2,
  polygon is a tuple or list of points"""
  a = 0
  for i in range( len( polygon)-1):
    a += angle_between_lines( point+polygon[i+1], point+polygon[i])
  a += angle_between_lines( point+polygon[0], point+polygon[ len( polygon)-1])
  a = abs( a)
  if a < 0.1:
    return False
  else:
    return True

def tkspline_to_quadratic_bezier( points):   #points = ((x1,y1),(x2,y2),...)
    if len(points) > 2:
      qbeziers = []
      startx,starty = points[0]             #start point
      curvex,curvey = points[1]             #point that defines curvature
      for nextx,nexty in points[2:-1]:
        endx = (curvex+nextx)/2             #calculate midway point
        endy = (curvey+nexty)/2             #to find end point
        qbeziers.append((startx,starty,curvex,curvey,endx,endy))
        
        startx,starty = endx,endy           #start where the last ended
        curvex,curvey = nextx,nexty
      
      endx,endy = points[-1]                #no midway calculation for last point
      qbeziers.append((startx,starty,curvex,curvey,endx,endy)) 
        
    else:
      startx,starty = points[0]             #If only two points are given
      endx,endy = points[1]                 #make a straight line
      curvex = (startx+endx)/2
      curvey = (starty+endy)/2
      qbeziers=((startx,starty,curvex,curvey,endx,endy),)
    
    return qbeziers

def tkspline_to_cubic_bezier( points):   #points = ((x1,y1),(x2,y2),...)
    qbeziers = tkspline_to_quadratic_bezier( points)
    cbeziers = []
    for qbez in qbeziers:
      startx, starty = qbez[0], qbez[1]
      endx, endy = qbez[4], qbez[5]
      qcurvex, qcurvey = qbez[2], qbez[3]
      
      ccurvexa, ccurveya = startx + 2/3*(qcurvex-startx), starty + 2/3*(qcurvey-starty)
      ccurvexb, ccurveyb = ccurvexa + 1/3*(endx-startx), ccurveya + 1/3*(endy-starty)
      cbeziers.append( (startx, starty, ccurvexa, ccurveya, ccurvexb, ccurveyb, endx, endy))
    # http://fontforge.sourceforge.net/bezier.html
    return cbeziers
    


def point_at_distance_from_line( x1, y1, x2, y2, d):
  "returns tuple of coordinates for a point in distance d from line orthogonal to point (x2,y2)"
  if round( y2, 3) -round( y1, 3) != 0:
    if y2 < y1:
      d = -d 
    k = -(x2-x1)/(y2-y1)
    x0 = ( d + sqrt( k**2 +1)*x2)/ sqrt( k**2 +1)
    y0 = y2 + k*( x0 -x2)
  else:
    if x1 > x2:
      d = -d
    x0 = x2
    y0 = y2 - d
  return (x0, y0)


def intersection_of_two_lines (x1,y1,x2,y2,x3,y3,x4,y4,parallel_detection_threshold=3): 
  """lines 1-2 and 3-4
  returns x,y, 0 if succesful or 1 if parallel
  y=mx+c is used
  parallel_detection_threshold is a negative decadic logarithm of minimal displacement of m
  that is considered parallel"""
  if x1-x2 == 0:                          
    if x3-x4 == 0:                      
      return 0,0,1,0        #lines paralell
    m2 = (y3-y4)/(x3-x4)
    c2 = y3 - m2 * x3
    rex,rey = x1,m2*x1+c2        
  elif x3-x4 == 0:
    m1 = (y1-y2)/(x1-x2)
    c1 = y1 - m1 * x1
    rex,rey = x3,m1*x3+c1
  else:
    m1 = (y1-y2)/(x1-x2)
    m2 = (y3-y4)/(x3-x4)
    if round(m1-m2,parallel_detection_threshold) == 0:
      return 0,0,1,0            #lines paralell
    c2 = y3 - m2 * x3
    c1 = y1 - m1 * x1
    rex = -(c2-c1)/(m2-m1)
    rey = (c1*m2-c2*m1)/(m2-m1)
  #check if point is on the lines
  online1 = is_point_beween_points_of_line( (x1,y1,x2,y2), (rex,rey))
  online2 = is_point_beween_points_of_line( (x3,y3,x4,y4), (rex,rey))
  if online1 and online2:
    online = 3
  elif online1:
    online = 1
  elif online2:
    online = 2
  else:
    online = 0
  # x-coord, y-coord , paralell(0 or 1), on line (0=no line, 1=on line 1-2, 2=on line 3-4, 3=on both) 
  return rex,rey,0,online


def mirror_point_on_line (xa,ya,x1,y1,x2,y2):
  """mirroring point a on line 1-2"""
  if x1-x2 == 0:                           #y=mc+c is used                                           
    return 2*x1-xa,ya
  elif y1-y2 == 0:
    return xa,2*y1-ya
  else:
    m=(y1-y2)/(x1-x2)
    cl=y1-m*x1
    ca=ya+(1/m)*xa
    xbl = (ca-cl)/(m+1/m)   #point bl is intersection with line
    ybl = -(1/m)*xbl+ca
    xb = xbl*2-xa
    yb = ybl*2-ya
    return xb,yb

def elongate_line (x1,y1,x2,y2,d):
  """line 1-2 will be elongatet at point 2 negative d will make it shorter"""
  if x1-x2 == 0:
    rex = x2
    if y2 > y1: 
      rey = y2 + d
    else:
      rey = y2 - d
  else:
    m = (y1-y2)/(x1-x2)
    dx = sqrt( d**2 / (1 + m**2))
    dy = m * dx
    if dy < 0:
      dy=-dy
    if d<0: 
      dx,dy = -dx, -dy
    if x2 > x1:
      rex = x2 + dx
    else:
      rex = x2 - dx
    if y2 > y1:
      rey = y2 + dy
    else:
      rey = y2 - dy
  return rex, rey


def coordinate_flat_list_to_xy_tuples( coords):
  ret = []
  for i in range( 0, len( coords), 2):
    ret.append( (coords[i],coords[i+1]))
  return ret
    

def find_parallel_polyline( coords, d):
  """take and returns list of tuples of coordinates for parallel abscissa in distance d"""
  res = []
  i = 0
  for x,y in coords:
    if i == 0:
      x1, y1 = coords[1]
      res.append( point_at_distance_from_line( x1,y1,x,y,-d))
    elif i == len( coords)-1:
      x1, y1 = coords[-2]
      res.append( point_at_distance_from_line( x1,y1,x,y,d))
    else:
      #internal point
      x0, y0 = coords[i-1] #previous point
      x1, y1 = coords[i+1] #next point
      xa1,ya1,xa2,ya2 = find_parallel( x0, y0, x, y, d)
      xb1,yb1,xb2,yb2 = find_parallel( x1, y1, x, y, -d)
      xx, yy, parallel, online = intersection_of_two_lines( xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2,parallel_detection_threshold=1)
      if not parallel:
        if online:
          # it the cross-point is on one of the lines, this could not be the case
          # where the point runs too far from the original line
          res.append( (xx, yy))
        else:
          # we should check for very near parallels going back -
          # the cross-point is very far from the original line, we need to fix it
          xa, ya = elongate_line(xa1,ya1,xa2,ya2,abs(d))
          if point_distance(xa,ya,xa2,ya2) < point_distance(xa,ya,xx,yy):
            xb, yb = elongate_line(xb1,yb1,xb2,yb2,abs(d))
            res.append( (xa, ya))          
            res.append( (xb, yb))
          else:
            res.append( (xx, yy))
      else:
        # is 1-2-3 straight or does 2-3 go back from 1-2
        if point_distance(x,y,x0,y0) > point_distance(x0,y0,x1,y1):
          d *= -1
          res.append((xa2,ya2)) # going back
        else:
          res.append((xa2,ya2)) # going forward
    i += 1
  return res



def point_on_quadratic_bezier( coords, t):
  """based on info from http://en.wikipedia.org/wiki/B%C3%A9zier_curve"""
  (x0,y0,x1,y1,x2,y2) = coords
  x = (1-t)**2*x0 + 2*t*(1-t)*x1 + t**2*x2
  y = (1-t)**2*y0 + 2*t*(1-t)*y1 + t**2*y2
  return x,y



def quadratic_beziere_to_polyline( point, n=10):
  t = 0
  points = []
  for i in range( n+1):
    points.append( point_on_quadratic_bezier( point, t))
    t += 1.0/n
  return points


def is_point_beween_points_of_line( line, point):
  """this function asserts that the point is on a line and
  computes if it is between the points defining the line"""
  x1, y1, x2, y2 = line
  x, y = point
  if point_distance( x1,y1,x,y) + point_distance( x2,y2,x,y) > 1.02 * point_distance( x1,y1,x2,y2):
    return False
  return True


def plane_normal_from_3_points( point1, point2, point3):
  for point in (point1,point2,point3):
    if None in point:
      return None  # some coords are missing
  from . import transform
  x1,y1,z1 = point1
  x2,y2,z2 = point2
  x3,y3,z3 = point3
  m0 = transform.matrix( [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])
  m1 = transform.matrix( [[1,y1,z1],[1,y2,z2],[1,y3,z3]])
  m2 = transform.matrix( [[x1,1,z1],[x2,1,z2],[x3,1,z3]])
  m3 = transform.matrix( [[x1,y1,1],[x2,y2,1],[x3,y3,1]])
  #d0 = 1.0*m0.get_determinant()
  d1 = 1.0*m1.get_determinant()
  d2 = 1.0*m2.get_determinant()
  d3 = 1.0*m3.get_determinant()
  a = d1 #/d0
  b = d2 #/d0
  c = d3 #/d0
  return a,b,c

def angle_between_planes( plane1, plane2):
  a1,b1,c1 = plane1
  a2,b2,c2 = plane2
  cos = (a1*a2 + b1*b2 + c1*c2) / sqrt( a1**2+b1**2+c1**2) / sqrt( a2**2+b2**2+c2**2)
  return cos

def same_or_oposite_side( plane1,plane2):
  a1,b1,c1 = plane1
  a2,b2,c2 = plane2
  a = a1+a2
  b = b1+b2
  c = c1+c2
  if (a**2 + b**2 + c**2) > (a1**2 + b1**2 + c1**2):
    return 1
  else:
    return -1


def expand_rectangle( coords, d):
  x1,y1,x2,y2 = coords
  return x1-d,y1-d,x2+d,y2+d


def create_transformation_to_coincide_point_with_z_axis( mov, point):
  """takes 3d coordinates 'point' (vector mov->point) and returns 3d transform object
  (transform3d.transform3d) that performs rotation to get 'point' onto z axis (x,y)=(0,0)
  with positive 'z'.
  NOTE: this is probably far from efficient, but it works
  """
  from .transform3d import transform3d
  t = transform3d()
  a,b,c = mov
  t.set_move( -a, -b, -c)
  x,y,z = t.transform_xyz( *point)
  t.set_rotation_y( atan2( x, z))
  x,y,z = t.transform_xyz( *point)
  t.set_rotation_x( -atan2( y, sqrt(x**2+z**2)))
  x,y,z = t.transform_xyz( *point)
  if z < 0:
      t.set_rotation_x( pi)
  #t.set_move( *mov)
  return t


def create_transformation_to_rotate_around_particular_axis( line_start, line_end, theta):
  """
  taken from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
  """
  a,b,c = line_start
  u,v,w = line_end
  u -= a
  v -= b
  w -= c
  u2 = u*u;
  v2 = v*v;
  w2 = w*w;
  cosT = cos(theta);
  sinT = sin(theta);
  l2 = u2 + v2 + w2;
  l =  sqrt(l2);

  if (l2 < 0.000000001):
      raise ValueError("RotationMatrix: direction vector too short!")

  m11 = (u2 + (v2 + w2) * cosT)/l2;
  m12 = (u*v * (1 - cosT) - w*l*sinT)/l2;
  m13 = (u*w * (1 - cosT) + v*l*sinT)/l2;
  m14 = (a*(v2 + w2) - u*(b*v + c*w) 
      + (u*(b*v + c*w) - a*(v2 + w2))*cosT + (b*w - c*v)*l*sinT)/l2;

  m21 = (u*v * (1 - cosT) + w*l*sinT)/l2;
  m22 = (v2 + (u2 + w2) * cosT)/l2;
  m23 = (v*w * (1 - cosT) - u*l*sinT)/l2;
  m24 = (b*(u2 + w2) - v*(a*u + c*w) 
      + (v*(a*u + c*w) - b*(u2 + w2))*cosT + (c*u - a*w)*l*sinT)/l2;

  m31 = (u*w * (1 - cosT) - v*l*sinT)/l2;
  m32 = (v*w * (1 - cosT) + u*l*sinT)/l2;
  m33 = (w2 + (u2 + v2) * cosT)/l2;
  m34 = (c*(u2 + v2) - w*(a*u + b*v) 
      + (w*(a*u + b*v) - c*(u2 + v2))*cosT + (a*v - b*u)*l*sinT)/l2;
  from .transform3d import transform3d
  t = transform3d( [[m11,m12,m13,m14],[m21,m22,m23,m24],[m31,m32,m33,m34],[0,0,0,1]])
  return t


  
  
