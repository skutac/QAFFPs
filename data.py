from __future__ import print_function
import os, csv

import numpy as np
import chembl_wrapper, utils

import config as cfg


def export_target_sets(db=None, DATA=None):
    chembl = chembl_wrapper.ChEMBLWrapper()
    targets = chembl.get_chembl_targets()
    count = len(targets)

    for i, t in enumerate(targets):
        print("\n{}/{}".format(i, count))
        bioactives = chembl.get_bioactives_for_target(t["tid"])
        bioactives = chembl.filter_by_activity_comment(bioactives)
        target_sets = chembl.export_compounds_for_target_by_activity_type(t, bioactives, 50)

def export_fingerprints_for_target_sets(TARGET_SETS=cfg.DIRS["TARGET_SETS"], fingerprints_dir=cfg.DIRS["FPS"]):
    target_sets = set(os.listdir(TARGET_SETS))
    done_fingerprints = set(os.listdir(fingerprints_dir))

    target_sets = list(target_sets - done_fingerprints)
    print("Export fingerprints for target sets: {}".format(len(target_sets)))

    for ts in target_sets:
        fps = []
        filepath = os.path.join(TARGET_SETS, ts)

        with open(filepath, "r") as input_file:
            mols = [m for m in csv.DictReader(input_file)]

        for m in mols:
            try:                
                fps.append((m["cmpd_chembl_id"], m["value"], utils.get_morgan_fingerprint_for_smiles(m["canonical_smiles"], radius=2, length=1024)))
            except Exception, e:
                print(e)

        export_fingerprints_to_file(os.path.join(fingerprints_dir, ts), fps)
    return

def export_fingerprints_to_file(filepath, fps):
    with open(filepath, "w") as output:
        csv_writer = csv.writer(output)
        csv_writer.writerow(["id", "value", "fp"])
        csv_writer.writerows(fps)
    return

def get_qsar_models_stats(indir=cfg.DIRS["QSAR_MODELS"]):
    files = [x for x in os.listdir(indir) if "csv" in x]

    chembl = chembl_wrapper.ChEMBLWrapper()
    tid2targets = {x["tid"]: x for x in chembl.get_chembl_targets()}
    result = []

    for f in files:
        filepath = os.path.join(indir, f)
        
        with open(filepath, "r") as r:
            rows = [x for x in csv.DictReader(r)]

        type2data = utils.group_records_by_field(rows, "set")
        
        tid = int(f.split("_")[0])
        row = {"target_id": tid, "activity_type": f.split(".")[0].split("_")[-1].lower(), "organism": tid2targets[tid]["organism"], "target_chembl_id": tid2targets[tid]["chembl_id"], "target_name": tid2targets[tid]["pref_name"]}
        
        for t in ["cv", "test"]:
            for p in ["rmse", "r20", "q2"]:
                values = [float(x[p]) for x in type2data[t]]
                label = "{}_{}".format(t if t == "test" else "train", p)
                row[label] = round(np.mean(values), 2)

        filepath = os.path.join(cfg.DIRS["TARGET_SETS"], f)
        
        with open(filepath, "r") as r:
            row["ligand_count"] = len([x for x in csv.DictReader(r)])

        result.append(row)


    header = ["target_id", "target_chembl_id", "target_name", "organism", "activity_type", "ligand_count", "train_rmse", "train_q2", "train_r20", "test_rmse", "test_q2", "test_r20"]
    
    with open(os.path.join(cfg.DIRS["RESULTS"], "qsar_models_stats.csv"), "w") as w:
        writer = csv.DictWriter(w, fieldnames=header)
        writer.writeheader()
        writer.writerows(result)
