from __future__ import print_function
import csv, os, math

import numpy

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.DataStructs import cDataStructs

import utils, config as cfg

def get_bqaffps(ligand_set_name, cutoff=5, confidence=90, max_dev=2, get_rdkit_obj=False):
    print("Generating b-QAFFPs for set {} (confidence: {}, max_dev: {})".format(ligand_set_name, confidence, max_dev))
    qaffps = get_qaffps(ligand_set_name, confidence=confidence, max_dev=max_dev)
    bqaffps = {}
    length = len(qaffps.values()[0])
    on_counts = []

    for cid, qaffp in qaffps.items():            
        bits = [1 if v is not None and v >= cutoff else 0 for v in qaffp]
        on_counts.append(bits.count(1))

        if get_rdkit_obj:
            fp = cDataStructs.ExplicitBitVect(length)
            on_indexes = [i for i, b in enumerate(bits) if b == 1]
            fp.SetBitsFromList(on_indexes)
            bits = fp
        
        bqaffps[cid] = bits

    print("Average of ON bits: {}".format(numpy.average(on_counts)))
    return bqaffps

def get_qaffps(ligand_set_name, confidence=90, max_dev=2, get_experimental=False, as_dicts=False):
    print("Get QAFFPs for {}".format(ligand_set_name))
    filepath = os.path.join(cfg.DIRS["QAFFPS"], "{}.csv".format(ligand_set_name))
    qaffps = {}

    if os.path.exists(filepath): 
        with open(filepath, "r") as input_file:
            original_rows = [r for r in csv.reader(input_file, delimiter=",")]
        
        header = original_rows[0]

        for r in original_rows[1:]:
            qaffps[r[0]] = [float(v) if v != "" else None for v in r[1:]]

    else:
        qaffps, header, stats = generate_qaffps(ligand_set_name, confidence=confidence, max_dev=max_dev, get_experimental=get_experimental)

    if as_dicts:
        dicts = []
        
        for r in rows[1:]:
            dicts.append({val:r[i]  for i, val in enumerate(rows[0])})

        rows = dicts
    
    return qaffps

def generate_qaffps(ligand_set_name, confidence=90, max_dev=1.0, get_experimental=False, export=True):
    print("Generating QAFFPS for set {} (confidence: {}, dev: {}, experimental: {})".format(ligand_set_name, confidence, max_dev, get_experimental))
    
    predictions = os.path.join(cfg.DIRS["PREDICTIONS"], ligand_set_name)
    prediction_files = os.listdir(predictions)
    prediction_files.sort(key = lambda x: (int(x.split("_")[0]), x))

    with open(os.path.join(predictions, prediction_files[0]), "r") as input_file:
        rows = [r for r in csv.DictReader(input_file)]        
        compound_set = list({c["id"] for c in rows})

    header = ["id"]
    compound2fingerprint = {}
    stats = {"compounds": len(compound_set), "confidence": confidence, "max_dev": max_dev, "original": 0, "none": 0, "predicted": 0}

    for f in prediction_files:
        target_set_id = f.split(".")[0]
        header.append(target_set_id)

        if get_experimental:
            with open(os.path.join(cfg.DIRS["TARGET_SETS"], f), "r") as input_file:
                compound2value = {r["cmpd_chembl_id"]:float(r["value"]) for r in csv.DictReader(input_file)}
        
        with open(os.path.join(predictions, f), "r") as input_file:
            compound2predicted = {r["id"]:r for r in csv.DictReader(input_file)}

        
        if confidence is not None:
            confidence, max_dev = float(confidence), float(max_dev)

            with open(os.path.join(cfg.DIRS["ERROR_MODELS"], f), "r") as error_values_file:
                alpha_values = [float(r["alpha"]) for r in csv.DictReader(error_values_file)]
                alpha_values.sort()
                pos = int(math.ceil(len(alpha_values)*(confidence/100.0))-1)
                alpha_value = alpha_values[pos]

        for c in compound_set:
            
            if get_experimental and c in compound2value:
                value = round(compound2value[c], 3)
                stats["original"] += 1
            else:
                if confidence is None:
                    value = round(float(compound2predicted[c]["value"]), 3)
                    stats["predicted"] += 1

                else:
                    dev = float(compound2predicted[c]["error"])*alpha_value

                    if dev <= max_dev:
                        value = round(float(compound2predicted[c]["value"]), 3)
                        stats["predicted"] += 1
                    else:
                        value = None
                        stats["none"] += 1
        
            if c in compound2fingerprint:
                compound2fingerprint[c].append(value)
            else:
                compound2fingerprint[c] = [value]

    print("Original values: {}".format(stats["original"]))
    print("Predicted values: {}".format(stats["predicted"]))
    print("None values: {}".format(stats["none"]))

    if export:
        if confidence is None:      
            outdir = "all"
        else:
            outdir = "{}_{}".format(confidence, max_dev)

        outpath = os.path.join(cfg.DIRS["QAFFPS"], outdir)
        
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        
        with open(os.path.join(outpath, "{}.csv".format(ligand_set_name)), "w") as output_file:
            writer = csv.writer(output_file)
            writer.writerow(header)

            for c, vals in compound2fingerprint.items():
                row = [c]
                row.extend(vals)
                writer.writerow(row)
    
    return compound2fingerprint, header, stats