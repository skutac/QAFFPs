import os

# ChEMBL database connection settings
DB = "chembl"
HOST = ""
USER = ""
PASSWD = ""

# RF parameters
N_ESTIMATORS = 100 # No. trees in RF
N_JOBS = 2 # No. used cores

# Directory settings
DATA_DIR = "data"
DIRS = {
    "DATA" : DATA_DIR,
    "LIGAND_SETS": os.path.join(DATA_DIR, "ligand_sets/"),
    "ERROR_MODELS" : os.path.join(DATA_DIR, "error_models/"),
    "FPS": os.path.join(DATA_DIR, "fingerprints/"),
    "QSAR_SETS": os.path.join(DATA_DIR, "qsar_sets/"),
    "QSAR_MODELS": os.path.join(DATA_DIR, "qsar_models/"),
    "RESULTS": os.path.join(DATA_DIR, "results/"),
    "TARGET_SETS": os.path.join(DATA_DIR, "target_sets/"),
    "PREDICTIONS": os.path.join(DATA_DIR, "predictions/"),
    "QAFFPS": os.path.join(DATA_DIR, "qaffps/"),
}

if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)

for d in DIRS.values():
    if not os.path.exists(d):
        os.makedirs(d)