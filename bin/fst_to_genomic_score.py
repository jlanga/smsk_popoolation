#!/usr/bin/env python3
import sys
from numpy import sqrt

pop1 = int(sys.argv[1])
pop2 = int(sys.argv[2])

def triangular_number(T):
    '''Finds if a number is triangular, and if so, computer the base of that triangle
    Necessary to know the number of populations involved in the analysis of the Fst'''
    n = ( 1 + sqrt( 1 + 8 * T ))/2
    if n.is_integer():
        return n
    print("The number of analysis made seem wrong. The number of F_st is not a triangular number. Exiting")
    sys.exit()

if pop1 == pop2:
    print("Populations must be different")
    sys.exit()

if pop1 >= pop2:
    pop1 , pop2 = pop2 , pop1 # Swap


for line in sys.stdin:
    line        = line.strip().split("\t")
    n_analysis  = len(line) - 5
    npop        = int(triangular_number( n_analysis ))
    chromosome  = line[0]
    position    = line[1]
    analysis    = line[5:]
    if pop2 >= npop:
        print("Error with analysis to be extracted, exitting")
        sys.exit()
    index = pop1 + pop2 - 1
    fst = analysis[ index ]
    fst = fst.split("=")[1]
    sys.stdout.write(chromosome + "\t" + position + "\t" + str(fst) + "\n")
