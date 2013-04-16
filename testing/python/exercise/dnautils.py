COMPLEMENTS = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def complement(sequence):
    """
    Calculate the complementary sequence of a DNA sequence.

    @param sequence: DNA sequence expressed as a lower-case string.
    @return complementary sequence.
    """
    cdna = ''
    try:
        for ch in sequence:
            cdna += COMPLEMENTS[ch]
        return cdna
    except TypeError:
        raise ValueError('The input is not a sequence e.g. a string or list')
    except KeyError:
        raise ValueError('The input is not a sequence of G,T,C,A')

def inverse(sequence):
    """
    Calculate the inverse of a DNA sequence.

    @param sequence: a DNA sequence expressed as an upper-case string.
    @return inverse as an upper-case string. 
    """
    # Reverse string using approach recommended on StackOverflow
    # http://stackoverflow.com/questions/931092/reverse-a-string-in-python
    return sequence[::-1]

def antiparallel(sequence):
    """
    Calculate the antiparallel of a DNA sequence.

    @param sequence: a DNA sequence expressed as an upper-case string.
    @return antiparallel as an upper-case string. 
    """
    return inverse(complement(sequence))
