import xgboost as xgb


def loadModel(model_path):
    model = xgb.XGBRanker()
    model.load_model(model_path)
    return model
