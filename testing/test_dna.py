from dna import calculate_weight

from dna import NUCLEOTIDES

def test_a():
    assert calculate_weight('A') == NUCLEOTIDES['A']

def test_g():
    assert calculate_weight('G') == NUCLEOTIDES['G']

def test_ga():
    assert calculate_weight('GA') == NUCLEOTIDES['G'] + NUCLEOTIDES['A']
