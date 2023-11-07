import os
import sys

from support.readFiles import readFiles
from support.train import train
from support.load_model import loadModel
from support.getacc import pred
from support.disableSomeFeature import transfer

# used_lst = [6, 7, 8, 9]
# all_num = [i for i in range(20)]
# lst = [i for i in all_num if i not in used_lst]
# print(lst)
# tmp = [".test", ".train", ".vali"]
# for i in tmp:
#     transfer("lp_3" + i, lst)
# exit(0)

if __name__ == "__main__":
    data_folder = "/home/yzz/ML/NewTraining_Tools_old/data"
    name = "lp_3"
    # get data
    all_data = readFiles(os.path.join(data_folder, name))

    # train model
    # train(all_data)

    # load model
    model_path = "/home/yzz/ML/NewTraining_Tools_old/model_1/cvrp_model_1.bin"
    model = loadModel(model_path)

    # get accuracy¬
    #  model, take_out, percent):
    take_out = [10, 15, 20]
    percent = [1, 0.9, 0.8, 0.7, 0.6]
    acc = pred(all_data, model, take_out, percent)
    for i in take_out:
        for k in range(len(percent)):
            print("take_out: {}, percent: {}, acc: {}".format(i, percent[k], acc[(i, k)]))
