import os


def create_model_path(base_path="model"):
    model_id = 0
    model_path = f"{base_path}_{model_id}"

    while os.path.exists(model_path):
        model_id += 1
        model_path = f"{base_path}_{model_id}"

    os.mkdir(model_path)
    return model_path
