import os.path

import pandas as pd
from sklearn.datasets import load_svmlight_file
from sklearn.datasets import dump_svmlight_file


def transfer(orig_data_path, removed_feature_set):
    # copy the data from the original file to a temporary file
    root = "data/"
    with open(os.path.join(root, orig_data_path), "r") as orig_file:
        tmp_name = "temp_{}".format(orig_data_path)
        tmp_path = os.path.join(root, tmp_name)
        with open(tmp_path, "w") as temp_file:
            for line in orig_file:
                temp_file.write(line)
    X, y = load_svmlight_file(tmp_path)
    data = pd.DataFrame(X.toarray())
    # print(data.iloc[:, 64:67])
    print(data.shape)
    # print(data)
    # always keep 0 and the last column (inferred from the original data)
    # safe_cols = [0, data.shape[1] - 1]
    safe_cols = [i for i in range(0, data.shape[1])]
    # print(safe_cols)
    data = data.drop(columns=removed_feature_set)
    # print(data.iloc[:, 65:67])
    libsvm_data = []
    df = data
    target = y
    col = list(df.columns)
    # if for every entry in safe_cols, the corresponding column is not in df.columns,
    # then add the column to df, but with all 0s
    # if all([i not in col for i in safe_cols]):
    for i in safe_cols:
        if i not in col:
            # add a column with all 0s to df
            df[i] = [0] * df.shape[0]
    # sort the columns by the column name
    # print(df)
    df = df.sort_index(axis=1)
    col = list(df.columns)
    # print(col)
    # print(df)
    # generate the libsvm data
    for index, row in df.iterrows():
        libsvm_line = f"{target[index]} "
        for i, value in enumerate(row):
            # if value != 0 or col[i] in safe_cols:
            libsvm_line += f"{col[i]}:{value} "
        libsvm_data.append(libsvm_line)
    all_data = "\n".join(libsvm_data)
    with open(tmp_path, "w") as temp_file:
        temp_file.write(all_data)

    # X = df.values
    # y = target
    # libsvm_data = pd.DataFrame({'target': y})
    # for i in range(X.shape[1]):
    #     nonzero = X[:, i] != 0
    #     libsvm_data[f'feature{i}'] = X[:, i]
    #     libsvm_data.loc[nonzero, f'feature{i}'] = f'{i}:{libsvm_data.loc[nonzero, f"feature{i}"].values}'
    # libsvm_data.to_csv('libsvm_data.txt', sep = ' ', header = False, index = False)
    #
    # # alternatively, you can use the dump_svmlight_file function from scikit-learn to write the libsvm data
    # dump_svmlight_file(X, y, tmp_path)


def addFeature(orig_data_path, add_feature_set):
    # copy the data from the original file to a temporary file
    root = "data/"
    with open(os.path.join(root, orig_data_path), "r") as orig_file:
        tmp_name = "temp_{}".format(orig_data_path)
        tmp_path = os.path.join(root, tmp_name)
        with open(tmp_path, "w") as temp_file:
            for line in orig_file:
                temp_file.write(line)
    X, y = load_svmlight_file(tmp_path)
    data = pd.DataFrame(X.toarray())
    # print(data.iloc[:, 64:67])
    print(data.shape)
    # print(data)
    # always keep 0 and the last column (inferred from the original data)
    safe_cols = [0, data.shape[1]]
    # print(safe_cols)

    # take the product of the columns in add_feature_set in one_more_col
    one_more_col = data[add_feature_set[0]] * data[add_feature_set[1]] * 1e6
    # add one more zero column to the data
    data[data.shape[1]] = one_more_col
    libsvm_data = []
    df = data
    target = y
    col = list(df.columns)
    # print(col)
    # if for every entry in safe_cols, the corresponding column is not in df.columns,
    # then add the column to df, but with all 0s
    if all([i not in col for i in safe_cols]):
        for i in safe_cols:
            if i not in col:
                # add a column with all 0s to df
                df[i] = [0] * df.shape[0]
    # sort the columns by the column name
    # print(df)
    df = df.sort_index(axis=1)
    col = list(df.columns)
    # print(df)
    # generate the libsvm data
    for index, row in df.iterrows():
        libsvm_line = f"{target[index]} "
        for i, value in enumerate(row):
            libsvm_line += f"{col[i]}:{value} "
        libsvm_data.append(libsvm_line)
    all_data = "\n".join(libsvm_data)
    with open(tmp_path, "w") as temp_file:
        temp_file.write(all_data)

    # X = df.values
    # y = target
    # libsvm_data = pd.DataFrame({'target': y})
    # for i in range(X.shape[1]):
    #     nonzero = X[:, i] != 0
    #     libsvm_data[f'feature{i}'] = X[:, i]
    #     libsvm_data.loc[nonzero, f'feature{i}'] = f'{i}:{libsvm_data.loc[nonzero, f"feature{i}"].values}'
    # libsvm_data.to_csv('libsvm_data.txt', sep = ' ', header = False, index = False)
    #
    # # alternatively, you can use the dump_svmlight_file function from scikit-learn to write the libsvm data
    # dump_svmlight_file(X, y, tmp_path)
