import os
import sys
import re
import json
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


BASE_PATH = '/data/archive/downstream/'
SOURCES = ['TARGET', 'TCGA', 'THR', 'TH']

AML = 'acute myeloid leukemia'
ALL = 'acute lymphoblastic leukemia'


def load_tsv(path):
    return pd.read_csv(path, sep='\t', index_col=0)

def donor_to_samples(donor_id):
    sample_paths = glob.glob(BASE_PATH + donor_id + '*')
    sample_ids = [p[p.rfind('/')+1:] for p in sample_paths]
    return sample_ids

def sample_to_donor(sample_id):
    src = sample_to_source(sample_id)
    if src == 'TARGET':
        return sample_id[:16]
    elif src == 'TCGA':
        return sample_id[:12]
    elif src == 'THR':
        return sample_id[:10]
    elif src == 'TH':
        return sample_id[:9]
    else:
        raise Exception('Unknown data source')

def sample_to_source(sample_id):
    for src in sorted(SOURCES, key=lambda s: -len(s)):
        if sample_id[:len(src)] == src:
            return src
    raise Exception('Unknown source for sample_id %s' % sample_id)

def sample_id_info(sample_id):
    info = {}
    for src in sorted(SOURCES, key=lambda s: -len(s)):
        if sample_id[:len(src)] == src:
            info['source'] = src
            break
    donor_lens = {'TARGET': 16, 'TCGA': 12, 'THR': 10, 'TH': 9, 'UNK': 0}
    info['donor'] = sample_id[donor_lens[info.get('source', 'UNK')]]
    m = re.search(r'[0-9][0-9]', sample_id)
    if m:
        info['sub_source'] = sample_id[0:m.end()]
    elif 'source' in info:
        info['sub_source'] = info['source']
    return info
    
def get_samples_by_donor_list(donor_ids):
    sample_ids = []
    missing = 0
    for donor_id in donor_ids:
        donor_samples = donor_to_samples(donor_id)
        if len(donor_samples):
            sample_ids.append(donor_samples[0])
        else:
            missing += 1
    print('Number of donors missing from file structure:', missing)
    return sample_ids

def get_samples_by_disease(disease, type_df, meta_df):
    donor_ids = list(type_df.loc[type_df['Diagnosis/Disease'] == disease].index)
    samples_from_type_df = get_samples_by_donor_list(donor_ids)
    samples_from_meta_df = list(meta_df.loc[meta_df['disease'] == disease].index)
    # there is some cross-over between the two
    return list(set(samples_from_type_df + samples_from_meta_df))


def get_sample_correlations(sample_id, corr_df):
    corr_dict = None
    
    # try to get correlations directly from all v all matrix
    if sample_id in corr_df.index:
        corr_row = corr_df.loc[sample_id]
        corr_dict = {s: c for s, c in zip(corr_df.index, corr_row)}
    
    # otherwise retrieve from file structure
    else:
        json_path = BASE_PATH + sample_id + '/tertiary/treehouse-8.0.1_comp-v5/2.0.json'
        if os.path.isfile(json_path):
            with open(json_path, 'r') as f:
                corr_dict = json.loads(f.read())['correlations_vs_focus_sample']
    
    return corr_dict

def get_samples_correlated_above_threshold(sample_id, corr_df, threshold=0.87):
    corr_dict = get_sample_correlations(sample_id, corr_df)
    # many samples do not have correlations either in the all v all matrix or file structure
    if corr_dict is None:
        return None
    
    sample_ids = []
    for s, c in corr_dict.items():
        if c > threshold:
            sample_ids.append(s)
    
    return sample_ids

def diseases_correlated_above_threshold(sample_id, corr_df, type_df, threshold=0.87, same_src=True, diff_src=True):
    sample_ids = get_samples_correlated_above_threshold(sample_id, corr_df, threshold)
    if sample_ids is None:
        return None
    
    result_sample_ids = []
    src = sample_to_source(sample_id)
    if same_src:
        result_sample_ids.extend([s for s in sample_ids if sample_to_source(s) == src])
    if diff_src:
        result_sample_ids.extend([s for s in sample_ids if sample_to_source(s) != src])
    sample_ids = result_sample_ids
    
    donor_ids = [sample_to_donor(sample_id) for sample_id in sample_ids]
    disease_counts = type_df.loc[donor_ids]['Diagnosis/Disease'].value_counts()
    return disease_counts / sum(disease_counts)
