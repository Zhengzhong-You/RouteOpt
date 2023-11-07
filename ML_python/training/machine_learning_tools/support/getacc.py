import heapq
import numpy as np


def pred(all_data, model, take_out, percent):
    _, _, _, _, \
        _, _, _, _, \
        data, _, sb_scores, group, = all_data
    length = len(group)
    percent_len = len(percent)
    beg = 0
    end = 0
    acc = {}
    cnt = {}
    for j in range(length):
        end += group[j]
        tmp_data = data[beg:end]
        pre_dict_y = np.array(model.predict(tmp_data))
        arr_large1 = [heapq.nlargest(min(take_out[i], len(pre_dict_y)), range(len(pre_dict_y)), pre_dict_y.take) for
                      i in
                      range(len(take_out))]
        arr_large = {}
        for i in range(len(take_out)):
            arr_large[take_out[i]] = arr_large1[i]
        sb = np.array(sb_scores[beg:end])
        max_score = heapq.nlargest(1, sb)[0]

        max_array = {}
        for i in take_out:
            val_ = [sb[j] for j in arr_large[i]]
            max_array[i] = heapq.nlargest(1, val_)[0]

        std = [int(i * max_score) for i in percent]
        for i in take_out:
            for k in range(percent_len):
                if max_array[i] >= std[k]:
                    acc[(i, k)] = acc.get((i, k), 0) + 1
                cnt[(i, k)] = cnt.get((i, k), 0) + 1
        beg = end

    for i in take_out:
        for k in range(percent_len):
            acc[(i, k)] /= cnt[(i, k)]

    return acc
