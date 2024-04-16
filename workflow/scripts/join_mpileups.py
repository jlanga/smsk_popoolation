#!/usr/bin/env python

import sys
import gzip
from functools import reduce


if __name__ == "__main__":
    files = sys.argv[1:]
    file_handlers = [gzip.open(f, "rt") for f in files]
    for lines in zip(*file_handlers):
        lines = [line.strip().split("\t") for line in lines]
        row_data = lines[0][0:3]
        sequence_data = reduce(lambda x, y: x + y, [line[3:6] for line in lines])
        sys.stdout.write("\t".join(row_data + sequence_data))
        sys.stdout.write("\n")
