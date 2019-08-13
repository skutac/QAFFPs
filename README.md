# QAFFPs

Files and scripts used for the generation of QAFFPs

First, you need to setup the connection to ChEMBL database (https://www.ebi.ac.uk/chembl/) running on a PostgreSQL server in config.py file, and setup the Python 2.7 virtual environment using the requirements.txt file (pip install -r requirements.txt).

You also need to install RDkit (https://www.rdkit.org/) with Python bindings, and the caret package for R (http://caret.r-forge.r-project.org/).

To run the workflow for the generation of QAFFPs, execute python run.py command in your console (see the code below).


```python
from __future__ import print_function
import os, subprocess
from subprocess import CalledProcessError, check_output

import chembl_wrapper, data, predictions, utils, config as cfg, qaffps as qaffps_lib

def run():
    # Export target sets from ChEMBL
    data.export_target_sets()

    # Export Morgan2 fingerprint for target sets
    data.export_fingerprints_for_target_sets()
    
    # Generate train/test sets     
    try:
        output = check_output(["Rscript", "--vanilla", "data_split.r", cfg.DIRS["FPS"], cfg.DIRS["QSAR_SETS"]])
        returncode = 0
    
    except CalledProcessError as e:
        output = e.output
        returncode = e.returncode
        raise Exception("ERROR IN data_split.r. See the output above.")

    # Build QSAR models
    target_sets = list({x.split(".")[0] for x in os.listdir(cfg.DIRS["TARGET_SETS"])})
    count = len(target_sets)

    for i, target_set in enumerate(target_sets, 1):
        print("{}/{}".format(i, count))
        predictions.build_qsar_models(target_set)

    # Get QSAR models stats
    data.get_qsar_models_stats()

    # Export Morgan2 fingerprints for a ligand set
    utils.prepare_ligand_set_from_set_file(os.path.join(cfg.DIRS["LIGAND_SETS"], "example_set.csv"))

    # Predict ligand set on all models
    predictions.predict_ligands_on_all_models(os.path.join(cfg.DIRS["FPS"], "example_set.csv"), r20_cutoff=0.6, q2_cutoff=0.5)
    
    # Generate QAFFPs for the ligand set
    qaffps_lib.generate_qaffps("example_set", confidence=90, max_dev=2)

    # Get QAFFPs
    qaffps = qaffps_lib.get_qaffps("example_set")

    #Get b-QAFFPs
    bqaffps = qaffps_lib.get_bqaffps("example_set", cutoff=5)
    
if __name__ == '__main__':
    run()
```