# Given two rectangles return a rectangle representing where
# they overlap or None if they do not overlap.
# Each rectangle is represented as a list [x0, x1, y0, y1] where
#  0 < x0 < x1
#  0 < y0 < y1

def overlap_axis(a0, a1, b0, b1):
    start = max(a0, b0)
    end = min(a1, b1)
    if (end <= start):
        return None
    return [start, end]

def overlap(a, b):
    a_x0, a_x1, a_y0, a_y1 = a
    b_x0, b_x1, b_y0, b_y1 = b
    overlap_x = overlap_axis(a_x0, a_x1, b_x0, b_x1)
    overlap_y = overlap_axis(a_y0, a_y1, b_y0, b_y1)
    if (overlap_x == None) or (overlap_y == None):
        return None
    return overlap_x + overlap_y
