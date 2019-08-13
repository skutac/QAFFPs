import os

HOST = ""
USER = "chembl"
PASSWD = "chembl"
DB = "chembl_25"

N_ESTIMATORS = 100
N_JOBS = 2

DATA = "data"
DIRS = {
    "DATA" : DATA,
    "LIGAND_SETS": os.path.join(DATA, "ligand_sets/"),
    "ERROR_MODELS" : os.path.join(DATA, "error_models/"),
    "FPS": os.path.join(DATA, "fingerprints/"),
    "QSAR_SETS": os.path.join(DATA, "qsar_sets/"),
    "QSAR_MODELS": os.path.join(DATA, "qsar_models/"),
    "RESULTS": os.path.join(DATA, "results/"),
    "TARGET_SETS": os.path.join(DATA, "target_sets/"),
    "PREDICTIONS": os.path.join(DATA, "predictions/"),
    "QAFFPS": os.path.join(DATA, "qaffps/"),
}

if not os.path.exists(DATA):
    os.makedirs(DATA)

for d in DIRS.values():
    if not os.path.exists(d):
        os.makedirs(d)