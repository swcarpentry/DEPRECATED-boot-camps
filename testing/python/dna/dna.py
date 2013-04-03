
NUCLEOTIDES = {'A':131.2, 'T':304.2, 'C':289.2, 'G':329.2}

def calculate_weight(sequence):
    """
    Calculate the molecular weight of a DNA sequence.

    @param sequence: DNA sequence expressed as an upper-case string.
    @return molecular weight.
    """
    weight = 0.0
    for ch in sequence:
        weight += NUCLEOTIDES[ch]
    return weight
