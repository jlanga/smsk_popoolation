#!/usr/bin/env python3
from numpy import mean
import sys


def get_blocks(stdin):
    """
    Retrieve records from stdin. Each record starts with a ">".
    """
    block = []
    for line in stdin:
        if block and line.startswith('>'):
            yield block
            block = []
        block.append(line.strip())
    if block:
        yield block


def compute_hp(block):
    """
    Compute the H_p
    """
    header = block.pop(0).split(" ")[0]
    header = header[1:]  # Trim the >
    header = header.split(":")
    chromosome, position = header[0:2]
    hps = []
    for line in block:

        if not line:
            break

        line = line.strip().split("\t")
        genotypes = map(int, line[4:8])
        genotypes = sorted(genotypes, reverse=True)
        maj_n, min_n = genotypes[0:2]
        hp = 2 * maj_n * min_n / ((maj_n + min_n)**2)
        hps.append(hp)

    result_string = "{chromosome}\t{position}\t{hp}\n".format(
        chromosome=chromosome,
        position=position,
        hp=mean(hps)
    )

    return result_string


if __name__ == "__main__":
    hps = (compute_hp(block) for block in get_blocks(sys.stdin))
    for result in hps:
        sys.stdout.write(result)
