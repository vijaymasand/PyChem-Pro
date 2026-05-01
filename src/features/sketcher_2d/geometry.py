# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from math import pi, atan2, cos, sin, asin, sqrt, degrees

# ---------------------------- POINT -------------------------------
def point_distance(p1, p2):
    """ calculate distance between two points """
    return sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)


def points_within_range(p1, p2, range_):
    """ check if p1 and p2 has x and y difference not higher than range_ """
    return abs(p1[0] - p2[0]) <= range_ and abs(p1[1] - p2[1]) <= range_


# ------------------------------- LINE --------------------------------
def line_contains_point(line, point):
    """ checks if line contains point """
    x1, y1, x2, y2 = line
    x, y = point
    if point_distance((x1, y1), (x, y)) + point_distance((x2, y2), (x, y)) > 1.02 * point_distance((x1, y1), (x2, y2)):
        return False
    return True


def line_length(line):
    return sqrt((line[2] - line[0])**2 + (line[3] - line[1])**2)


def line_extend_by(line, d):
    """ returns the point upto which it gets extended.
    -ve value of d will make it shorter """
    x1, y1, x2, y2 = line
    if x1 - x2 == 0:
        rex = x2
        if y2 > y1:
            rey = y2 + d
        else:
            rey = y2 - d
    else:
        m = (y1 - y2) / (x1 - x2)
        dx = sqrt(d**2 / (1 + m**2))
        dy = m * dx
        if dy < 0:
            dy = -dy
        if d < 0:
            dx, dy = -dx, -dy
        if x2 > x1:
            rex = x2 + dx
        else:
            rex = x2 - dx
        if y2 > y1:
            rey = y2 + dy
        else:
            rey = y2 - dy
    return rex, rey


def line_get_point_at_distance(line, d):
    """ returns a point at perpendicular distance d from line's end point.
    -ve value of d gives point from opposite side (right side on screen) """
    x1, y1, x2, y2 = line
    if round(y2, 3) - round(y1, 3) != 0:
        if y2 < y1:
            d = -d
        k = -(x2 - x1) / (y2 - y1)
        x0 = (d + sqrt(k**2 + 1) * x2) / sqrt(k**2 + 1)
        y0 = y2 + k * (x0 - x2)
    else:
        if x1 > x2:
            d = -d
        x0 = x2
        y0 = y2 - d
    return x0, y0


def line_get_parallel(line, d):
    """ returns the parallel line at distance d.
    -ve value of d gives line from opposide side (right side on screen) """
    x1, y1, x2, y2 = line
    if round(y2, 3) - round(y1, 3) != 0:
        if y2 < y1:
            d = -d
        k = -(x2 - x1) / (y2 - y1)
        x = (d + sqrt(k**2 + 1) * x1) / sqrt(k**2 + 1)
        y = y1 + k * (x - x1)
        x0 = (d + sqrt(k**2 + 1) * x2) / sqrt(k**2 + 1)
        y0 = y2 + k * (x0 - x2)
    else:
        if x1 > x2:
            d = -d
        x, x0 = x1, x2
        y = y1 - d
        y0 = y2 - d
    return x, y, x0, y0


def line_get_intersection_of_line(line, line2):
    """ returns x,y, status. (status -> 0=successful, 1=parallel) """
    x1, y1, x2, y2 = line
    x3, y3, x4, y4 = line2
    parallel_detection_threshold = 3

    if x1 - x2 == 0:
        if x3 - x4 == 0:
            return 0, 0, 1  # parallel
        m2 = (y3 - y4) / (x3 - x4)
        c2 = y3 - m2 * x3
        rex, rey = x1, m2 * x1 + c2
    elif x3 - x4 == 0:
        m1 = (y1 - y2) / (x1 - x2)
        c1 = y1 - m1 * x1
        rex, rey = x3, m1 * x3 + c1
    else:
        m1 = (y1 - y2) / (x1 - x2)
        m2 = (y3 - y4) / (x3 - x4)
        if round(m1 - m2, parallel_detection_threshold) == 0:
            return 0, 0, 1  # parallel
        c2 = y3 - m2 * x3
        c1 = y1 - m1 * x1
        rex = -(c2 - c1) / (m2 - m1)
        rey = (c1 * m2 - c2 * m1) / (m2 - m1)
    return rex, rey, 0


def line_get_side_of_point(line, point, threshold=0):
    """ tells whether a point is on one side of a line or on the other.
    return vals are 1=left, -1=right, 0=point on line. (sides are for screen coordinate) """
    x1, y1, x2, y2 = line
    x, y = point
    a = atan2(y - y1, x - x1)
    b = atan2(y2 - y1, x2 - x1)
    if a * b < 0 and abs(a - b) > pi:
        if a < 0:
            a += 2 * pi
        else:
            b += 2 * pi
    if abs(a - b) <= threshold or abs(abs(a - b) - pi) <= threshold:
        return 0
    elif a - b < 0:
        return 1
    else:
        return -1


def line_get_angle_from_east(line):
    """ returns the angle between the center-east line and 'line'.
    angle is clockwise on screen """
    angle = atan2(line[3] - line[1], line[2] - line[0])
    if angle < 0:
        angle += 2 * pi
    return angle


# ---------------------------- RECTANGLE ---------------------------
def rect_get_center(rect):
    return (rect[0] + rect[2]) / 2, (rect[1] + rect[3]) / 2


def rect_contains_point(rect, pt):
    x1, y1, x2, y2 = rect_normalize(rect)
    return x1 <= pt[0] <= x2 and y1 <= pt[1] <= y2


def rect_normalize(rect):
    x1, y1, x2, y2 = rect
    if x2 < x1:
        x2, x1 = x1, x2
    if y2 < y1:
        y2, y1 = y1, y2
    return [x1, y1, x2, y2]


def rect_intersects_rect(rect1, rect2):
    xs = rect1[0], rect1[2], rect2[0], rect2[2]
    ys = rect1[1], rect1[3], rect2[1], rect2[3]
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    w1 = abs(rect1[0] - rect1[2])
    h1 = abs(rect1[1] - rect1[3])
    w2 = abs(rect2[0] - rect2[2])
    h2 = abs(rect2[1] - rect2[3])
    if w1 + w2 > dx and h1 + h2 > dy:
        return True
    return False


def rect_get_intersection_of_line(rect, line):
    """ returns the intersection of line and rectangle's perimeter.
    if no intersection found, returns None. """
    x1, y1, x2, y2 = rect
    lx1, ly1, lx2, ly2 = line
    # Lines forming the rectangle
    rect_lines = [
        (x1, y1, x2, y1),
        (x2, y1, x2, y2),
        (x2, y2, x1, y2),
        (x1, y2, x1, y1)
    ]
    intersections = []
    for r_line in rect_lines:
        res_x, res_y, status = line_get_intersection_of_line(line, r_line)
        if status == 0:
            # Check if intersection point is within both line segments
            # We use a small epsilon for robustness
            if line_contains_point(r_line, (res_x, res_y)):
                # Also check distance from line start to avoid returning the center
                intersections.append((res_x, res_y))
    
    if not intersections:
        return None
        
    # Pick the intersection closest to (lx1, ly1) but still on the perimeter
    # (Usually there's only one if lx1 is inside, or two if lx1 is outside)
    intersections.sort(key=lambda p: (p[0] - lx1)**2 + (p[1] - ly1)**2)
    return intersections[0]


# ---------------------------- CIRCLE --------------------------------
def circle_get_point(center, radius, direction, resolution=0):
    dx, dy = direction[0] - center[0], direction[1] - center[1]
    if resolution:
        angle = round(atan2(dy, dx) / (pi * resolution / 180.0)) * (pi * resolution / 180.0)
    else:
        angle = atan2(dy, dx)
    x = center[0] + round(cos(angle) * radius, 2)
    y = center[1] + round(sin(angle) * radius, 2)
    return x, y


def rotate_point(px, py, cx, cy, angle):
    dx = px - cx
    dy = py - cy
    rot_x = dx * cos(angle) - dy * sin(angle)
    rot_y = dx * sin(angle) + dy * cos(angle)
    return cx + rot_x, cy + rot_y


# Matrix Multiplication
def matrix_multiply_3(A, B):
    AB = []
    for i in range(len(A)):
        AB.append([])
        for j in range(len(B[0])):
            AB[i].append(A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j])
    return AB


class Transform:
    def __init__(self, mat=None):
        if mat:
            self.mat = mat
            return
        self.mat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    def transform(self, x, y):
        x_res, y_res, w = matrix_multiply_3(self.mat, [[x], [y], [1]])
        return x_res[0], y_res[0]

    def scale(self, scale):
        self.mat = matrix_multiply_3([[scale, 0, 0], [0, scale, 0], [0, 0, 1]], self.mat)

    def rotate(self, angle):
        self.mat = matrix_multiply_3([[cos(angle), -sin(angle), 0], [sin(angle), cos(angle), 0], [0, 0, 1]], self.mat)

    def translate(self, tx, ty):
        self.mat = matrix_multiply_3([[1, 0, tx], [0, 1, ty], [0, 0, 1]], self.mat)
