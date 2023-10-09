import os
import sys
import re


def generateIdxIns(path):
    idxIns = []
    path1 = "./../../../../" + path
    files = os.listdir(path1)
    for file in files:
        if file.endswith(".vrp"):
            idxIns.append("./../../../" + path + file)
    with open("OldIns.ins", 'w') as f:
        for file in idxIns:
            f.write(file + "\n")


generateIdxIns("DataForCVRP/OldIns/")
