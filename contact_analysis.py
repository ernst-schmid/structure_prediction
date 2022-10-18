#!/usr/bin/python

from argparse import ArgumentParser
from importlib import machinery
import pdb_parser3 as pdb_parser
import json, sys, math, re
import statistics as stat

import shutil
import os
from os import listdir,mkdir
from os.path import isfile, join

def get_complete_complexes(path, str = '.done.txt'):

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    complex_names = []
    for name in onlyfiles:
        if str in name:
            complex_names.append(name[:-9])
    
    return complex_names

def get_files_for_complex(path, complex_name, filetype = ''):

    return [f for f in listdir(path) if complex_name in f and f.endswith(filetype)]

def process_pdb_data(path, pdb_filenames, angstrom_cutoff, plddt_cutoff, pae_cutoff):

    for pdb_filename in pdb_filenames:

        model_num = re.findall(r'model_\d+', pdb_filename)[0]
        scores_filename = pdb_filename.replace('.pdb', '_scores.json')
        scores_file = open(join(path, scores_filename), 'r')
        score_data = scores_file.read()
        score_data = json.loads(score_data)
        pae_data = score_data['pae']
        plddt_data = score_data['plddt']
        scores_file.close()

        contacts = pdb_parser.get_contacts(join(path, pdb_filename), angstrom_cutoff)
        final_contacts = []
        for c in contacts:

            aa1 = c['aa1'] 
            aa2 = c['aa2']
            ix1 = aa1['abs_index']
            ix2 = aa2['abs_index']

            pae_val = pae_data[ix1-1][ix2-1]
            aa1['plddt'] = plddt_data[ix1-1]
            aa2['plddt'] = plddt_data[ix2-1]
            if(pae_val > pae_cutoff or aa1['plddt'] < plddt_cutoff or aa2['plddt'] < plddt_cutoff):
                continue

            contact_id = aa1['type'] + str(ix1) + aa2['type'] + str(ix2)
            
            if contact_id in final_contacts:
                final_contacts[contact_id]['instances'].append({"data":c, "pae":pae_val, "found_in_model":model_num})
            else:
                final_contacts[contact_id]= {'instances':[{"data":c, "pae":pae_val, "found_in_model":model_num}]}

    return final_contacts
        


def process_folder(folder_path, angstrom_cutoff):

    complex_names = get_complete_complexes(folder_path)
    print("Found " + str(len(complex_names)) + " files to process");

    num_complete = 0
    for complex_name in complex_names:

        pdb_filenames = get_files_for_complex(folder_path, complex_name, '.pdb')
        process_pdb_data(folder_path, pdb_filenames, angstrom_cutoff)

        num_complete += 1
        print("Finished (" + str(round(100*num_complete/len(complex_names), 1)) + "%) " + complex_name);
