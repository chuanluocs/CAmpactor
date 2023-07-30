"""
This file is used for translating the output of AutoCCAG, FastCA and TCA
into the format of output of SamplingCA. 
"""

import os, sys

def format_translate(from_path, to_path):
    lis = []

    with open(from_path, "r") as f:
        for l in f.readlines():
            line = l.strip()
            if line:
                curlis = [str((int(x) & 1)) for x in line.split()]
                lis.append(curlis)

    with open(to_path, "w") as f:
        for curlis in lis:
            f.write(" ".join(curlis) + "\n")
