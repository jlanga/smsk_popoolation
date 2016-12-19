#!/usr/bin/env python
from numpy import mean
import sys


def get_records(stdin):
    block = []
    for line in stdin:
        if block and line.startswith('>'):
            yield block
            block = []
        block.append(line.strip())
    if block:
        yield block

# def mount_header(line):
#     line = line.strip()


def compute_Hp(lines):
    header     = lines.pop(0).split(" ")[0]
    header     = header[1:] # Trim the >
    header     = header.split(":")
    chromosome = header[0]
    position   = header[1]
    HP = []
    for line in lines:
        
        if not line:
            break
        
        line = line.strip().split("\t")
        
        genotypes  = map( int , line[4:8] )
        
        genotypes = sorted( genotypes, reverse = True)
        maj_n   = genotypes[0]
        min_n   = genotypes[1]
        hp    = 2 * maj_n * min_n / ( ( maj_n + min_n )^2 )
        HP.append(hp)
        
    return chromosome + "\t" + position + "\t" + str(mean(HP)) + "\n"

for record in get_records(sys.stdin):
    sys.stdout.write( compute_Hp( record ) )
