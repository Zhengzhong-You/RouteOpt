from support import convertRank
import random

split = [0.5, 0.3, 0.2]


# float(p.split(':')[1]) != 0.0 and

def save_data(group_data, all_data, mapping, root_tag):
    train_output_feature, train_output_SB_scores, train_output_group, vali_output_feature, vali_output_SB_scores, \
    vali_output_group, test_output_feature, test_output_SB_scores, test_output_group, \
    root_br_feature, root_br_SB_scores, root_br_group = all_data

    if len(group_data) == 0:
        return
    SB_scores = [float(data[0]) for data in group_data]
    if_discard = convertRank.convertRank(group_data, mapping)
    if if_discard:
        return
    a = random.uniform(0, 1)
    if a < split[0]:
        train_output_group.write(str(len(group_data)) + "\n")
        cnt = 0
        for data in group_data:
            # only include nonzero features
            feats = [p for p in data[2:]]
            train_output_feature.write(data[0] + " " + " ".join(feats) + "\n")
            train_output_SB_scores.write(str(SB_scores[cnt]) + "\n")
            cnt += 1
    elif a < (split[0] + split[1]):
        vali_output_group.write(str(len(group_data)) + "\n")
        cnt = 0
        for data in group_data:
            # only include nonzero features
            feats = [p for p in data[2:]]
            vali_output_feature.write(data[0] + " " + " ".join(feats) + "\n")
            vali_output_SB_scores.write(str(SB_scores[cnt]) + "\n")
            cnt += 1
    else:
        test_output_group.write(str(len(group_data)) + "\n")
        cnt = 0
        for data in group_data:
            # only include nonzero features
            feats = [p for p in data[2:]]
            test_output_feature.write(data[0] + " " + " ".join(feats) + "\n")
            test_output_SB_scores.write(str(SB_scores[cnt]) + "\n")
            cnt += 1
        if root_tag:
            root_br_group.write(str(len(group_data)) + "\n")
            cnt = 0
            for data in group_data:
                # should include zero features in root_br
                feats = [p for p in data[2:]]
                root_br_feature.write(data[0] + " " + " ".join(feats) + "\n")
                root_br_SB_scores.write(str(SB_scores[cnt]) + "\n")
                cnt += 1
