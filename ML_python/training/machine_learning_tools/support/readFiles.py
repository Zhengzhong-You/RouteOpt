from sklearn.datasets import load_svmlight_file
import os


def readFiles(name):
    all_name = ["{}.train", "{}.vali", "{}.test"]
    files = []
    for i in all_name:
        transfer = i.format(name)
        x, y = load_svmlight_file(transfer)
        sb_scores = []
        with open("{}.SB_scores".format(transfer), "r") as f:
            data = f.readlines()
            for line in data:
                sb_scores.append(float(line.split("\n")[0]))
        group = []
        with open("{}.group".format(transfer), "r") as f:
            data = f.readlines()
            for line in data:
                group.append(int(line.split("\n")[0]))
        files.append(x)
        files.append(y)
        files.append(sb_scores)
        files.append(group)

    return files
