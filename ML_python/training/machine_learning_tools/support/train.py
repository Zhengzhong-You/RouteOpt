import xgboost as xgb
from support.getModelPath import create_model_path
from support.getacc import pred
import os


def train(all_data):
    x_train, y_train, _, group_train, \
        x_valid, y_valid, _, group_valid, \
        _, _, _, _ = all_data

    params = {
        'objective': 'rank:pairwise',
        'n_estimators': 1000,
        'tree_method': 'gpu_hist',
        'gpu_id': '0'
    }
    model = xgb.XGBRanker(**params)
    model.fit(x_train, y_train, group_train, verbose=False, eval_set=[(x_train, y_train), (x_valid, y_valid)],
              eval_group=[group_train, group_valid])

    model_path = create_model_path()
    model.save_model('{}/xgb.bin'.format(model_path))

    take_out = [10, 15, 20]
    percent = [1, 0.9, 0.8, 0.7, 0.6]
    acc = pred(all_data, model, take_out, percent)
    for i in take_out:
        for k in range(len(percent)):
            print("take_out: {}, percent: {}, acc: {}".format(i, percent[k], acc[(i, k)]))
    with open('{}/model_info.txt'.format(model_path), 'w') as f:
        f.write(str(params))
        f.write('\n')
        for i in take_out:
            for k in range(len(percent)):
                f.write("take_out: {}, percent: {}, acc: {}".format(i, percent[k], acc[(i, k)]))
                f.write('\n')
