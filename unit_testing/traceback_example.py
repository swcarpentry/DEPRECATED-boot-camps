import math
import sys

def area_of_a_circle(radius):
    '''
    Calculate the area of a circle
    '''
    return 2.0 * math.pi * radius

def area_of_a_square(half_side):
    '''
    Calculate the area of a square
    '''
    return (2.0*half_side) **2

def diff_area_circle_area_square(radius):
    '''
    Calcualte how much bigger the area of a square is than
    the area of a circle where the radius of the circle is
    half the length of one side of the square.
    '''
    area_circ = area_of_a_circle(radius)
    area_sq = area_of_a_square(radius)

    return area_sq - area_circ

if __name__ == "__main__":
    r = sys.argv[1]
    diff_area_circle_area_square(r)
