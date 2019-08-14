from __future__ import print_function
import csv, os, cPickle

import numpy

from sklearn import model_selection as cv, metrics
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

import validation, config as cfg, utils

def get_k_stratified_sets(input_set, k, y_index=2):
    input_set_size = len(input_set)
    input_set.sort(key = lambda x: x[y_index])
    folds = [[] for x in xrange(k)]

    index = 0
    
    for item in input_set:
        folds[index].append(item)
        index = index + 1
        if index%k == 0:
            index = 0

    return folds

def cross_validation(train, k):
    print("\nCross-validation: {}-fold".format(k))
    folds = get_k_stratified_sets(train, k)
    fold2data = {}

    for fold_index, fold in enumerate(folds):
        fold2data[fold_index] = {"y_true": [], "y_pred": [], "residuals": []}
        train_folds = [f for x, f in enumerate(folds) if x != fold_index]
        test_fold = fold

        train_set = []
        for f in train_folds:
            train_set.extend(f)

        train_x = [t[1] for t in train_set]
        train_y = [t[2] for t in train_set]
        test_x = [t[1] for t in test_fold]
        test_y = [t[2] for t in test_fold]

        rf = RandomForestRegressor(n_estimators=cfg.N_ESTIMATORS, n_jobs=cfg.N_JOBS)
        rf.fit(train_x, train_y)
        pred_y = rf.predict(test_x)

        fold2data[fold_index]["y_true"] = test_y
        fold2data[fold_index]["y_pred"] = pred_y
        fold2data[fold_index]["residuals"] = abs(pred_y-test_y)
        
        val_metrics = validation.ValidationMetrics(pred_y, test_y)        

        fold2data[fold_index]["rmse"] = val_metrics.RMSE
        fold2data[fold_index]["r20"] = val_metrics.Rsquared0
        fold2data[fold_index]["q2"] = val_metrics.Qsquared
        fold2data[fold_index]["ids"] = [t[0] for t in test_fold]

    for fold_index, fold in enumerate(folds):
        train_folds = [f for x, f in enumerate(folds) if x != fold_index]
        train_residuals = [f["residuals"] for x, f in fold2data.items() if x != fold_index]
        test_fold = fold

        train_set = []
        train_y = []
        for i, f in enumerate(train_folds):
            train_set.extend(f)
            train_y.extend(train_residuals[i])

        train_x = [t[1] for t in train_set]
        test_x = [t[1] for t in test_fold]
        test_y = fold2data[fold_index]["residuals"]

        cp = RandomForestRegressor(n_estimators=cfg.N_ESTIMATORS, n_jobs=cfg.N_JOBS)
        cp.fit(train_x, train_y)
        pred_y = cp.predict(test_x)
        alphas = abs(test_y/pred_y)
        fold2data[fold_index]["residuals_pred"] = pred_y
        fold2data[fold_index]["alphas"] = abs(alphas)

    for metric in ["RMSE", "R20", "q2"]:
        print("{}: {}".format(metric, numpy.mean([fold2data[i][metric.lower()] for i, fold in fold2data.items()])))

    return fold2data

def random_forest(train, test, repeat_count=10):
    print("RANDOM FOREST: {}(train), {}(test), {} repeats\n".format(len(train), len(test), repeat_count))
    train_x = [t[1] for t in train]
    train_y = [t[2] for t in train]
    test_x = [t[1] for t in test]
    test_y = [t[2] for t in test]
    max_r20 = -10
    repeat2data = {}

    for x in range(repeat_count):
        rf = RandomForestRegressor(n_estimators=cfg.N_ESTIMATORS, n_jobs=cfg.N_JOBS)
        rf.fit(train_x, train_y)

        if len(test):
            pred_y = rf.predict(test_x)
            
            val_metrics = validation.ValidationMetrics(pred_y, test_y)
            repeat2data[x] = {"rmse": val_metrics.RMSE, "r20": val_metrics.Rsquared0, "q2": val_metrics.Qsquared}
            repeat2data[x]["ids"] = [t[0] for t in test]
            repeat2data[x]["y_pred"] = pred_y

            if val_metrics.Rsquared0 > max_r20:
                max_r20 = val_metrics.Rsquared0
                final_rf = rf

        else:
            final_rf = rf

    return final_rf, repeat2data

def build_qsar_models(target_set_id):
    print("QSAR models for {}".format(target_set_id))
    train_ids = get_target_train_set(target_set_id)
    test_ids = get_target_test_set(target_set_id)
    fps = get_fps_for_target_set(target_set_id)
    print("Number of compounds: {}".format(len(fps)))

    train_set = [m for m in fps if m["id"] in train_ids]
    test_set = [m for m in fps if m["id"] in test_ids]
    train_list = [[m["id"], [int(bit) for bit in m["fp"]], m["value"]] for m in train_set]
    test_list = [[m["id"], [int(bit) for bit in m["fp"]], m["value"]] for m in test_set]

    fold2data = cross_validation(train_list, 10)
    
    print("\nCalculating QSAR model...")
    qsar_model, repeat2data = random_forest(train_list, test_list, 10)
    store_qsar_model(qsar_model, target_set_id, cfg.DIRS["QSAR_MODELS"])
    store_modeling_stats(fold2data, repeat2data, target_set_id)
    
    molid2residual = {}
    for i, fold in fold2data.items():
        for x, molid in enumerate(fold["ids"]):
            molid2residual[molid] = fold["residuals"][x]
    
    train_list_cp = []
    for item in train_list:
        train_list_cp.append([item[0], item[1], molid2residual[item[0]]])

    print("Calculating ERROR model...")
    conformal_model, repeat2data = random_forest(train_list_cp, [], 1)
    store_qsar_model(conformal_model, target_set_id, cfg.DIRS["ERROR_MODELS"])
    store_alphas(fold2data, target_set_id)
     
def get_target_test_set(target_set_id):
    with open(os.path.join(cfg.DIRS["QSAR_SETS"], "{}_test.csv".format(target_set_id)), "r") as input_file:
        rows = [r[0] for r in csv.reader(input_file)]
    return rows

def get_target_train_set(target_set_id):
    with open(os.path.join(cfg.DIRS["QSAR_SETS"], "{}_train.csv".format(target_set_id)), "r") as input_file:
        rows = [r[0] for r in csv.reader(input_file)]
    return rows

def get_fps_for_target_set(target_set_id):
    with open(os.path.join(cfg.DIRS["FPS"], "{}.csv".format(target_set_id)), "r") as input_file:
        rows = [{"id":r["id"], "fp":r["fp"], "value": float(r["value"])} for r in csv.DictReader(input_file)]
    return rows

def get_compound_set(target_set_id):
    with open(os.path.join(cfg.DIRS["TARGET_SETS"], "{}.csv".format(target_set_id))) as input_file:
        rows = [r for r in csv.DictReader(input_file)]
    return rows

def get_fingerprint_set(target_set_id):
    with open(os.path.join(cfg.DIRS["FPS"], "{}.csv".format(target_set_id))) as input_file:
        rows = [r for r in csv.DictReader(input_file)]
    return rows

def store_qsar_model(qsar_model, filename, outdir):
    model_file = open(os.path.join(outdir, "{}.cpickle".format(filename)), "w")
    cPickle.dump(qsar_model, model_file)

def store_alphas(fold2data, target_set_id):
    header = ["id", "res", "pred_res", "alpha"]
    output = csv.writer(open(os.path.join(cfg.DIRS["ERROR_MODELS"], "{}.csv".format(target_set_id)), "w"), delimiter=",")
    output.writerow(header)

    for i, fold in fold2data.items():
        for j, molid in enumerate(fold["ids"]):
            row = [molid, round(fold["residuals"][j], 3), round(fold["residuals_pred"][j], 3), round(fold["alphas"][j], 3)]
            output.writerow(row)

def store_modeling_stats(fold2data, repeat2data, filename):
    header = ["set", "rmse", "r20", "q2"]
    output = csv.writer(open(os.path.join(cfg.DIRS["QSAR_MODELS"], "{}.csv".format(filename)), "w"), delimiter=",")
    output.writerow(header)

    for i, fold in fold2data.items():
        row = []
        row.append("cv")
        row.extend([round(fold[h], 3) for h in header[1:]])
        output.writerow(row)

    for i, repeat in repeat2data.items():
        row = []
        row.append("test")
        row.extend([round(repeat[h], 3) for h in header[1:]])
        output.writerow(row)

def get_model(target_set_id):
    with open(os.path.join(cfg.DIRS["QSAR_MODELS"], "{}.cpickle".format(target_set_id))) as input_file:
        model = cPickle.load(input_file)
    return model

def get_error_model(target_set_id):
    with open(os.path.join(cfg.DIRS["ERROR_MODELS"], "{}.cpickle".format(target_set_id))) as input_file:
        model = cPickle.load(input_file)
    return model

def predict_ligands_on_all_models(ligands_file, r20_cutoff=0.6, q2_cutoff=0.5, force=True):
    with open(ligands_file, "r") as input_file:
        ligands = [r for r in csv.DictReader(input_file)]

    ligands_dir = ligands_file.split("/")[-1].split(".")[0]
    models = ["_".join([t["target_id"], t["activity_type"]]) for t in utils.get_models_by_parameters(r20_cutoff=r20_cutoff, q2_cutoff=q2_cutoff)]

    if not force:
        utils.check_dir(os.path.join(cfg.DIRS["PREDICTIONS_DIR"], ligands_dir))
        done = {x.split(".")[0] for x in os.listdir(os.path.join(cfg.DIRS["PREDICTIONS_DIR"], ligands_dir))}
        models = list(set(models).difference(done))
    
    count = len(models)
    print("Number of selected models: {}".format(count))

    fps = [[bit for bit in r["fp"]] for r in ligands]
    
    for m in models:
        pred_y = predict_fingerprints_on_model(fps, m)
        pred_error = predict_fingerprints_on_error_model(fps, m)
        data = [{"id": l["id"], "value": round(pred_y[i], 3), "error": round(pred_error[i], 3)} for i, l in enumerate(ligands)]
        utils.store_data_as_csv(os.path.join(cfg.DIRS["PREDICTIONS"], ligands_dir, "{}.csv".format(m)), data, ["id", "value", "error"])


def predict_fingerprints_on_model(fps, model_name):
    print("Predict {} fingerprints on model {}".format(len(fps), model_name))
    model = get_model(model_name)
    return model.predict(fps)

def predict_fingerprints_on_error_model(fps, model_name):
    print("Predict {} fingerprints on error model {}".format(len(fps), model_name))
    error_model = get_error_model(model_name)
    return error_model.predict(fps)
