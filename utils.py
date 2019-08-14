from __future__ import print_function
import csv, os

import psycopg2
from psycopg2 import extras

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.DataStructs import cDataStructs

import config as cfg

def postgres_array(items, mode="single"):
    if mode == "single":
        result = ",".join(["'{}'".format(i) for i in items])
    elif mode == "list":
        result = []
        for item in items:
            result.append(",".join(["'{}'".format(i) for i in item]))
        
        result = ",".join(["({})".format(i) for i in result])

    return result

def check_dir(path):
    if len(path.split(".")) > 1:
        path = "/".join(path.split("/")[:-1])
    
    if not os.path.exists(path):
        os.makedirs(path)

def group_records_by_field(records, field):
    field2records = {}
    
    if type(field) is list:
        
        for r in records:
            key = "_".join([str(r[k]) for k in field])
            if key in field2records:
                field2records[key].append(r)
            else:
                field2records[key] = [r]
        
    else:

        for r in records:
            key = r[field]
            if key in field2records:
                field2records[key].append(r)
            else:
                field2records[key] = [r]
    
    return field2records

def dictfetchall(cursor):
    """Return all rows from a cursor as a dict"""
    columns = [col[0] for col in cursor.description]
    return [
        dict(zip(columns, row))
        for row in cursor.fetchall()
    ]

def get_cursor():
    conn = psycopg2.connect(dbname=cfg.DB, user=cfg.USER, password=cfg.PASSWD, cursor_factory=psycopg2.extras.DictCursor)
    cur = conn.cursor()
    return cur

def get_morgan_fingerprint_obj_for_smiles(smiles, radius=2, length=1024):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length)
    return fp

def get_morgan_fingerprint_for_smiles(smiles, radius=2, length=1024):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length)
    fp = fp.ToBitString()
    return fp

def get_morgan_fingerprint_for_mol(mol, radius=2, length=1024):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length)
    return fp.ToBitString()

def get_models_by_parameters(r20_cutoff=0.6, q2_cutoff=0.5):
    with open(os.path.join(cfg.DIRS["RESULTS"], "qsar_models_stats.csv"), "r") as input_file:
        model_parameters = [r for r in csv.DictReader(input_file)]

    filtered_models = [m for m in model_parameters if float(m["train_q2"]) >= q2_cutoff and float(m["test_r20"]) >= r20_cutoff]    
    print("QSAR models [R20 >= {}, q2 >= {}]:".format(r20_cutoff, q2_cutoff), len(filtered_models))
    return list(filtered_models)

def prepare_ligand_set_from_set_file(input_file, output_file):
    with open(input_file, "r") as r:
        ligands = [row for row in csv.DictReader(r)]
    
    fps = []
    for l in ligands:
        mol = Chem.MolFromSmiles(l["smiles"])
        if mol:
            fp = get_morgan_fingerprint_for_mol(mol)
            fps.append([l["id"], fp])

    with open(output_file, "w") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(["id", "fp"])
        writer.writerows(fps)

    return fps

def store_data_as_csv(filename, datadicts, header=False):
    check_dir(filename)
    if not header:
        header = datadicts[0].keys()

    with open(filename, "w") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=header)
        writer.writeheader()
        writer.writerows(datadicts)
