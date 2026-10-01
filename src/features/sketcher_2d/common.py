# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)

def list_difference(list_):
    """return a list of differences between list members,
    the list is by 1 shorter than the original"""
    result = []
    for i in range(len(list_) - 1):
        result.append(list_[i] - list_[i + 1])
    return result


def difference(a, b):
    """ returns (a-b), i.e items of a which are not in b """
    ret = list(a)  # needed for type conversion of tuple for instance
    for i in b:
        if i in ret:
            ret.remove(i)
    return ret


def filter_unique(self, l):
    return list(dict.fromkeys(l))  # Requires python >= 3.7


def flatten(l):
    return [x for y in l for x in y]


def bbox_of_bboxes(bboxes):
    if len(bboxes) == 0:
        return
    xs = []
    ys = []
    for (x1, y1, x2, y2) in bboxes:
        xs += [x1, x2]
        ys += [y1, y2]
    return [min(xs), min(ys), max(xs), max(ys)]


# converts float to string by stripping trailing zeros.
# we can neither use simply ("%f" % num) as it creates trailing zeros,
# nor str(num) as it converts 0.00001 to undesired 1e-5
def float_to_str(num):
    # num is converted to float, so that it also works for int containing trailing zeros
    return str(float(round(num, 4))).rstrip("0").rstrip(".")


def str_to_tuple(strng):
    """ converts tuple str like '(1, 2, 3)' to tuple """
    return tuple(map(float, strng[1:-1].split(", ")))


def find_matching_parentheses(text, index):
    # @index is index of first parentheses
    # returns the index of last parentheses or index of last character
    # if no matching parentheses found
    size = len(text)
    i = index + 1
    count = 1
    while i < size and count:
        char = text[i]
        i += 1
        if char == "(":
            count += 1
        if char == ")":
            count -= 1
    return i - 1


def gen_combinations_of_series(series):
    """series is a list of lists (tuples), the generator yields
    lists by combining each element of each list with each other"""
    counter = len(series) * [0]
    end = [i - 1 for i in map(len, series)]
    counter[0] = -1  # dirty trick
    while counter != end:
        for i, e in enumerate(end):
            if counter[i] < e:
                counter[i] += 1
                for j in range(i):
                    counter[j] = 0
                break
        yield [s[counter[j]] for j, s in enumerate(series)]
