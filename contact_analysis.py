#!/usr/bin/python

from argparse import ArgumentParser
from importlib import machinery
import pdb_parser3 as pdb_parser
import orjson, json, sys, math, re
import statistics as stat
import numpy as np

import shutil
import os
from os import listdir,mkdir
from os.path import isfile, join

def get_complete_complexes(path, search_str = '.done.txt'):

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    complex_names = []
    for name in onlyfiles:
        if search_str in name:
            complex_names.append(name[:-9])
    
    return complex_names

def pdockq(x):

    return 0.724/(1 + 2.718**(-0.052*(x - 152.611))) + 0.018

def get_files_for_complex(path, complex_name, filetype = ''):

    return [f for f in listdir(path) if complex_name in f and f.endswith(filetype)]

def get_pae_data(scores_filename):
    scores_file = open(scores_filename, 'r')
    fstr = scores_file.read()
    pae_index = fstr.find('"pae":')
    pae_data = fstr[pae_index + 6:fstr.find(']]', pae_index) + 2].replace('[','').replace(']','').split(',')
    scores_file.close()

    return pae_data

def process_pdb_data(path, pdb_filenames, angstrom_cutoff, plddt_cutoff, pae_cutoff):

    aa_indices = {};
    final_contacts = []
    for pdb_filename in pdb_filenames:

        contacts = pdb_parser.get_contacts(join(path, pdb_filename), angstrom_cutoff)
        if len(contacts) < 3: 
            continue 
        model_num = re.findall(r'model_\d+', pdb_filename)[0]
        scores_filename = pdb_filename.replace('.pdb', '_scores.json')
        pae_data = None

        for c in contacts:
            
            aa1 = c['aa1'] 
            aa2 = c['aa2']
            ix1 = aa1['abs_index']
            ix2 = aa2['abs_index']

            if(aa1['plddt'] < plddt_cutoff or aa2['plddt'] < plddt_cutoff):
                continue

            if pae_data is None :
                pae_data = get_pae_data(join(path, scores_filename))

            pae_index = int(math.sqrt(len(pae_data)))*(ix1-1) + ix2-1
            pae_val = float(pae_data[pae_index])
            if(pae_val > pae_cutoff):
                continue

            chain_contact_id = aa1['chain'] + aa2['chain']
            if(chain_contact_id not in aa_indices):
                aa_indices[chain_contact_id] = [[], []]
            
            aa_indices[chain_contact_id][0].append(ix1)
            aa_indices[chain_contact_id][1].append(ix2)

            contact_id = aa1['chain'] + aa1['type'] + str(ix1) + aa2['chain'] + aa2['type'] + str(ix2)
            
            if contact_id in final_contacts:
                final_contacts[contact_id]['instances'].append({"data":c, "pae":pae_val, "found_in_model":model_num})
            else:
                final_contacts[contact_id]= {'instances':[{"data":c, "pae":pae_val, "found_in_model":model_num}]}

    return final_contacts

def process_folder(folder_path, angstrom_cutoff, plddt_cutoff, pae_cutoff, output_name):

    complex_names = get_complete_complexes(folder_path)
    print("Found " + str(len(complex_names)) + " files to process");

    num_complete = 0
    for complex_name in complex_names:

        pdb_filenames = get_files_for_complex(folder_path, complex_name, '.pdb')
        contacts = process_pdb_data(folder_path, pdb_filenames, angstrom_cutoff, plddt_cutoff, pae_cutoff)

        num_complete += 1
        print("Finished (" + str(round(100*num_complete/len(complex_names), 1)) + "%) " + complex_name);


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="output_structs",
    help="Directory with structure outputs, pdb_files, from Colabfold",
)

parser.add_argument(
    "--distance-cutoff",
    default=10,
    help="Minimum required distance in angstroms for a contact to be included in the final output ",
    type=float,
)

parser.add_argument(
    "--pae-cutoff",
    default=12,
    help="Minimum PAE value required for a contact to be included in the final output",
    type=float,
)

parser.add_argument(
    "--plddt-cutoff",
    default=70,
    help="Minimum plDDT values required by both residues in a contact in order for that contact to be included in the final output",
    type=float,
)

parser.add_argument(
    "--output-name",
    default="output",
    help="Name of the resulting output file in which to store the contacts",
)

args = parser.parse_args()

process_folder(args.input, args.distance_cutoff, args.plddt_cutoff, args.pae_cutoff, args.output_name)