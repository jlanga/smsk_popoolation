#!/usr/bin/env bash
# $1 = chromosome
# $2 = position of the center of the window
# $5 = tajimaD

awk '{ print $1"\t"$2"\t"$5 }' "$@"
