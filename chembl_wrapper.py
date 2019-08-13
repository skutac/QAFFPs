from __future__ import print_function

import csv, os, numpy, copy

import utils, config

ACCEPTED_TYPES = {
    "ic50":{"nM": {"unit": "-log(M)", "type": "pIC50"}},
    "pic50":{"None": {"unit": "-log(M)", "type": "pIC50"}},
    "log ic50":{"None": {"unit": "-log(M)", "type": "pIC50"}},
    "logic50":{"None": {"unit": "-log(M)", "type": "pIC50"}},

    "ec50":{"nM": {"unit": "-log(M)", "type": "pEC50"}},
    "pec50":{"None": {"unit": "-log(M)", "type": "pEC50"}},
    "log ec50":{"None": {"unit": "-log(M)", "type": "pEC50"}},
    "-log ec50":{"None": {"unit": "-log(M)", "type": "pEC50"}},

    "ki":{"nM": {"unit": "-log(M)", "type": "pKi"}},
    "pki":{"None": {"unit": "-log(M)", "type": "pKi"}},
    "log ki":{"None": {"unit": "-log(M)", "type": "pKi"}},
    "-log ki":{"none": {"unit": "-log(M)", "type": "pKi"}},

    "kd":{"nM": {"unit": "-log(M)", "type": "pKd"}},
    "log kd":{"nM": {"unit": "-log(M)", "type": "pKd"}},
    "pkd":{"None": {"unit": "-log(M)", "type": "pKd"}},
}

class ChEMBLWrapper():

    def __init__(self):
        self.cursor = utils.get_cursor()
        self.targets = False

    def get_chembl_targets(self, target_type="", keyword="", organism=""):
        """
        Get all targets from ChemblDB which match keyword and organism filters, if both None get them all
        """
        self.cursor.execute("""SELECT target_dictionary.tid, 
                                    target_type, 
                                    pref_name, 
                                    target_dictionary.chembl_id,
                                    target_dictionary.organism, 
                                    component_sequences.accession,
                                    protein_family_classification.l1,
                                    protein_family_classification.l2,
                                    protein_family_classification.l3
                            FROM target_dictionary
                            JOIN target_components ON target_components.tid = target_dictionary.tid
                            JOIN component_sequences ON component_sequences.component_id = target_components.component_id
                            JOIN component_class on component_class.component_id = component_sequences.component_id
                            JOIN protein_family_classification on protein_family_classification.protein_class_id = component_class.protein_class_id
                            WHERE target_type LIKE '%%%s%%'
                            AND pref_name LIKE '%%%s%%'
                            AND target_dictionary.organism LIKE '%%%s%%';""" %(target_type, keyword, organism))
        
        self.targets = self.preprocess_targets(utils.dictfetchall(self.cursor))
        print("Number of targets:", len(self.targets))
        return self.targets

    def preprocess_targets(self, targets):
        for i in xrange(len(targets)):
            if not targets[i]["organism"]:
                targets[i]["organism"] = "Unknown"

        return targets

    def get_chembl_target(self, tid):
        """
        Get Chembl target by its ID
        """
        self.cursor.execute("""SELECT * FROM target_dictionary
                            WHERE tid = '%s';""" %(tid))

        target = self.cursor.fetchall()
        if len(target) == 1:
            return target
        else:
            return False

    def get_bioactives_for_target(self, tid):
        """utils.postgres_array(std_cids)
        Get activities for given target
        """
        print("""Target ID:""", tid)
        self.cursor = utils.get_cursor()
        
        self.cursor.execute("""SELECT activities.activity_id,
            activities.pchembl_value,
            activities.assay_id,
            activities.molregno,
            activities.standard_relation ,
            activities.standard_value,
            activities.standard_units,
            activities.standard_flag,
            activities.standard_type,
            activities.activity_comment,
            assays.confidence_score,
            assays.assay_type,
            target_dictionary.target_type as tgt_type,
            target_dictionary.pref_name as tgt_pref_name,
            target_dictionary.chembl_id as tgt_chembl_id,
            target_dictionary.organism,
            molecule_dictionary.pref_name as cmpd_pref_name,
            molecule_dictionary.molecule_type,
            molecule_dictionary.structure_type,
            molecule_dictionary.inorganic_flag,
            compound_structures.canonical_smiles,
            compound_structures.standard_inchi_key,
            molecule_dictionary.chembl_id AS cmpd_chembl_id,
            assays.tid
            FROM activities,
            assays,
            target_dictionary,
            target_components,
            molecule_dictionary,
            compound_properties,
            compound_structures
            WHERE pchembl_value IS NOT NULL
            AND activities.data_validity_comment IS NULL
            AND activities.assay_id = assays.assay_id
            AND assays.tid = target_dictionary.tid
            AND target_dictionary.tid = target_components.tid
            AND molecule_dictionary.molregno = activities.molregno
            AND molecule_dictionary.molregno = compound_properties.molregno
            AND molecule_dictionary.molregno = compound_structures.molregno
            AND activities.standard_relation = '='
            AND assays.confidence_score IN ('7', '9')
            AND lower(activities.standard_type) IN ({})
            AND assays.tid = {}""".format(utils.postgres_array(ACCEPTED_TYPES.keys()), tid))

        compounds = utils.dictfetchall(self.cursor)
        print("Raw data:", len(compounds), "activities")
        return compounds

    def get_target_data(self, tid):
        """
        Get data for target by target_id
        """
        target = False
        if self.targets:
            for t in self.targets:
                if t["tid"] == tid:
                    target = t
                    break
        if not target:
            target = self.get_chembl_target(tid)

        return target

    def merge_activities(self, compounds):
        """Merge values for identical compounds - average of values, if STD > 0.5 reject the value"""

        chemblid2records = utils.group_records_by_field(compounds, "cmpd_chembl_id")
        merged = []
        
        for c, cs in chemblid2records.items():
            merged_record = chemblid2records[c][0]
            if len(cs) > 1:
                
                confidence_scores = list(set([r["confidence_score"] for r in cs]))
                if len(confidence_scores) > 1:
                    confidence_scores.sort()
                    highest_confidence = confidence_scores[-1]
                    merged_record["confidence_score"] = highest_confidence
                    values = [v["pchembl_value"] for v in chemblid2records[c] if v["confidence_score"] == highest_confidence]
                else:
                    values = [v["pchembl_value"] for v in chemblid2records[c]]
                
                values = [float(v) for v in values]
                values.sort()

                max_diff = 0.5
                sd = numpy.std(values)

                if sd <= max_diff: #the activities for one compound should be quite close together 
                    merged_record["value"] = numpy.average(values)
                    merged_record["merged"] = True
                    merged.append(merged_record)
                
            else:
                merged_record["merged"] = False
                merged_record["value"] = merged_record["pchembl_value"]
                merged.append(merged_record)

        print("Processed activities:", len(merged), "\n")
        return merged

    def filter_by_activity_comment(self, compounds):
        exclude = ["not active", "inconclusive", "inactive"]
        filtered = []
        
        for c in compounds:
            check = True

            if c["activity_comment"]:
                comment = c["activity_comment"].lower()
                for e in exclude:
                    if e in comment:
                        check = False
                        break
                if check:
                    filtered.append(c)
            else:
                filtered.append(c)
        return filtered


    def normalize_types(self, compounds):
        """Filter accepted activity types and convert them to standard units and types"""
        normalized = []
        for c in compounds:
            if c["standard_units"] is None:
                c["standard_units"] = "None"

            try:                
                standard_unit = c["standard_units"]
                standard_type = c["standard_type"].lower()

                # ONLY pChEMBL values accepted from ChEMBL 25
                # if standard_type in ["IC50", "Ki", "EC50", "Kd"] and c["standard_units"] == "nM" and not c["pchembl_value"]:
                #     value = round(-math.log(float(c["value"])/math.pow(10,9), 10), 2)
                #     c["value"] = value

                c["standard_type"] = ACCEPTED_TYPES[standard_type][standard_unit]["type"]
                c["standard_units"] = ACCEPTED_TYPES[standard_type][standard_unit]["unit"]
                normalized.append(c)

            except KeyError:
                pass

        print("After type normalization:", len(normalized))
        return normalized

    def export_data_set_to_file(self, compounds, target_id, activity_type):
        filepath = os.path.join(config.DIRS["TARGET_SETS"], "{}_{}.csv".format(target_id, activity_type.lower()))

        with open(filepath, "w") as output:
            csv_writer = csv.DictWriter(output, fieldnames = ["tid", "cmpd_chembl_id", "canonical_smiles", "standard_type", "value", "confidence_score", "cmpd_pref_name", "activity_comment"], extrasaction="ignore")
            csv_writer.writeheader()
            csv_writer.writerows(compounds)

    def export_compounds_for_target_by_activity_type(self, target, compounds, compound_treshold):
        """
        Initialization of data upload for one target by its activity types
        """

        compounds = copy.deepcopy(compounds)
        compounds = self.normalize_types(compounds)
        type2compounds = utils.group_records_by_field(compounds, "standard_type")
        result = {}

        for t, cs in type2compounds.items():
            if len(cs) >= compound_treshold:
            
                print("Activity type: {} \n{} compounds\n****************\n".format(t, len(cs)))
                cleaned_compounds = self.merge_activities(cs)

                if len(cleaned_compounds) > compound_treshold:
                    result[t] = cleaned_compounds
                    self.export_data_set_to_file(cleaned_compounds, target["tid"], t)

        return result