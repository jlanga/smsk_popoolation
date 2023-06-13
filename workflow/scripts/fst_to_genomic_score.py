#!/usr/bin/env python3

"""Script to extract the fst between two populations.

Usage: python3 fst_to_genomic_score.py 1 2 < popoolation2.fst > pop1_2.fst
Output is chr, pos, fst
"""

import sys
from numpy import sqrt

def triangular_number(number):
    '''Finds if a number is triangular, and if so, computer the base of that
    triangle. It is necessary to know the number of populations involved in the
    analysis of the Fst.

    >>> triangular_number(1)
    1
    >>> triangular_number(3)
    3
    >>> triangular_number(2)
    -1
    '''
    base = (1 + sqrt(1 + 8 * number))/2
    if base.is_integer():
        return int(base)
    return -1


def populations_to_index(pop1, pop2, npop):
    """Get the index of the column of a pair of populations in the fst file from
    popoolation2.

    Columns in the fst file are
    chromosome, center of the window, number of snps, covered fraction,
    coverage, fst 1:2, fst 1:3, ..., fst 1:n, fst 2:3, ..., fst 2:n, ...,
    fst  n-1:n

    Got from here: https://stackoverflow.com/a/27088560/4730336

    >>> populations_to_index(1, 2, 2)
    1
    >>> populations_to_index(1, 2, 22)
    1
    >>> populations_to_index(2, 3, 22)
    22
    >>> populations_to_index(21, 22, 22)
    231
    """
    pop1 = pop1 - 1
    pop2 = pop2 - 1
    big_triangle = int(npop * (npop - 1) / 2)
    small_triangle = int((npop - pop1) * (npop - pop1 - 1) / 2)
    k = big_triangle - small_triangle + pop2 - pop1
    return k


if __name__ == '__main__':

    POPULATION_1 = int(sys.argv[1])
    POPULATION_2 = int(sys.argv[2])

    if POPULATION_1 == POPULATION_2:
        raise SystemExit("Populations must be different")

    # Swap
    if POPULATION_1 >= POPULATION_2:
        POPULATION_1, POPULATION_2 = POPULATION_2, POPULATION_1

    NPOP = None

    for line in sys.stdin:
        line = line.strip().split("\t")
        n_analysis = len(line) - 5
        if not NPOP:
            NPOP = int(triangular_number(n_analysis))
        chromosome, position, _, _, _, *analysis = line
        if POPULATION_2 >= NPOP + 1:
            raise SystemExit("Error with analysis to be extracted, exitting")
        # index = pop1 + pop2 - 1
        index = populations_to_index(POPULATION_1, POPULATION_2, NPOP) - 1
        #sys.stderr.write(f"{index}")

        fst = analysis[index]
        fst = fst.split("=")[1]
        sys.stdout.write(f"{chromosome}\t{position}\t{fst}\n")
