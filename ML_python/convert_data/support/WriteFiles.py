import os
import support.SaveData as SaveData
import support.emptyFiles as emptyFiles


def writeFiles(tmp_dir, name, mapping):
    emptyFiles.cleanEmptyFile(tmp_dir)

    train_output_feature = open("data/{}_{}.train".format(name, mapping), "w")
    train_output_SB_scores = open("data/{}_{}.train.SB_scores".format(name, mapping), "w")
    train_output_group = open("data/{}_{}.train.group".format(name, mapping), "w")
    vali_output_feature = open("data/{}_{}.vali".format(name, mapping), "w")
    vali_output_SB_scores = open("data/{}_{}.vali.SB_scores".format(name, mapping), "w")
    vali_output_group = open("data/{}_{}.vali.group".format(name, mapping), "w")
    test_output_feature = open("data/{}_{}.test".format(name, mapping), "w")
    test_output_SB_scores = open("data/{}_{}.test.SB_scores".format(name, mapping), "w")
    test_output_group = open("data/{}_{}.test.group".format(name, mapping), "w")

    root_br_feature = open("data/{}_{}.root_br".format(name, mapping), "w")  # unseen data!
    root_br_SB_scores = open("data/{}_{}.root_br.SB_scores".format(name, mapping), "w")
    root_br_group = open("data/{}_{}.root_br.group".format(name, mapping), "w")

    all_data = [train_output_feature, train_output_SB_scores, train_output_group,
                vali_output_feature, vali_output_SB_scores, vali_output_group,
                test_output_feature, test_output_SB_scores, test_output_group,
                root_br_feature, root_br_SB_scores, root_br_group]

    files = os.listdir(tmp_dir)
    for file in files:
        if 'txt' in file:
            group_data = []
            group = ""
            with open(os.path.join(tmp_dir, file), 'r') as f:
                root_tag = False
                for line in f:
                    if not line:
                        break
                    if "#" in line:
                        line = line[:line.index("#")]
                    splits = line.strip().split(" ")
                    if splits[1] != group:
                        if len(group_data) != 0:
                            SaveData.save_data(group_data, all_data, mapping, root_tag)
                            root_tag = False
                        else:
                            root_tag = True
                        group_data = []
                    group = splits[1]
                    group_data.append(splits)
                SaveData.save_data(group_data, all_data, mapping, root_tag)
    for file in all_data:
        file.close()
