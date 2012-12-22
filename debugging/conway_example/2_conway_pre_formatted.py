#!/usr/bin/env python
"""
Conway's game of life example, part two.

This does not conform to PEP8 standards. Use pep8 to check most standards, and
then fix the others by eye.
"""


def conway(population, 
    generations = 100):
    for i in range(generations): population = evolve(population)
    return population

def evolve(population):
    activeCells = population | set(neighbor for p in population for neighbor in neighbors(p))
    newPopulation = set()
    for cell in activeCells:
        count = sum(neighbor in population for neighbor in neighbors(cell))
        if count == 3 or (count == 2 and cell in population): newPopulation.add(cell)
    return newPopulation

def neighbors(cell):
    x, y = cell
    return [(x, y), (x+1, y), (x-1, y), (x, y+1), (x, y-1), (x+1, y+1), (x+1, y-1), (x-1, y+1), (x-1, y-1)]

glider = set([(0, 0), (1, 0), (2, 0), (0, 1), (1, 2)])
print conway(glider)
