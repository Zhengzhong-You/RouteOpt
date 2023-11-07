import os


def cleanEmptyFile(path):
    tmp_dir = path
    files = os.listdir(tmp_dir)
    for file in files:
        if 'txt' in file:
            path = "{}/{}".format(tmp_dir, file)
            s = round(os.path.getsize(path) / float(1024), 2)
            if s < 0.01:
                os.remove(path)
