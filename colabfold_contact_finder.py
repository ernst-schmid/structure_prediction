#!/usr/bin/python

from argparse import ArgumentParser
from importlib import machinery
import pdb_parser
import json, sys, math, re
import statistics as stat

import shutil
import os
from os import listdir,mkdir
from os.path import isfile, join

def mean(values):

    if len(values) == 0:
        return 0
    sum = 0
    for d in values:
        sum += d

    return sum/len(values)


def pdoq_score(x):

    return 0.724/(1 + 2.718**(-0.052*(x - 152.611))) + 0.018
    

def get_all_files_for_complex(path, complex_name, filetype = ''):

    files_found = [f for f in listdir(path) if complex_name in f and f.endswith(filetype)]
    return files_found

def process_folder(path, angstrom_cutoff = 10, pae_cutoff = 14, plddt_cutoff = 70, output_filename = 'output', put_in_folders = False):

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    complex_names = []
    for name in onlyfiles:
        if '.done.txt' in name:
            complex_names.append(name[:-9])
    
    print("Found " + str(len(complex_names)) + " files to process");

    summary_output_data = []

    output_txt = 'complex_name\taa1\taa1_index\taa1_atom\taa1_plddt\taa2\taa2_index\taa2_atom\taa2_plddt\tangstrom_dist\tpae\tperc_models_found_in\tmodels_found_in\n'

    num_complete = 0
    for complex_name in complex_names:

        pdb_file_names = get_all_files_for_complex(path, complex_name, 'pdb')
        final_contacts = {};
        
        for pdb_filename in pdb_file_names:
            contacts = pdb_parser.get_contacts(join(path,pdb_filename), angstrom_cutoff)
            scores_filename = pdb_filename.replace('.pdb', '_scores.json')
            scores_file = open(join(path, scores_filename), 'r')
            score_data = scores_file.read()
            score_data = json.loads(score_data)
            pae_data = score_data['pae']
            plddt_data = score_data['plddt']
            scores_file.close()

            model_num = re.findall(r'model_\d+', pdb_filename)[0]

            for contact in contacts:

                aa1 = contact['aa1'] 
                aa2 = contact['aa2']
                ix1 = aa1['abs_index']
                ix2 = aa2['abs_index']

                pae_val = pae_data[ix1-1][ix2-1]
                aa1['plddt'] = plddt_data[ix1-1]
                aa2['plddt'] = plddt_data[ix2-1]
                if(pae_val > pae_cutoff or aa1['plddt'] < plddt_cutoff or aa2['plddt'] < plddt_cutoff):
                    continue

                contact_id = aa1['type'] + str(ix1) + aa2['type'] + str(ix2)
                
                if contact_id in final_contacts:
                    final_contacts[contact_id]['instances'].append({"data":contact, "pae":pae_val, "found_in_model":model_num})
                else:
                    final_contacts[contact_id]= {'instances':[{"data":contact, "pae":pae_val, "found_in_model":model_num}]}
        
        plddt_averages = [];
        aa1_indices = []

        for contact_id, contact in final_contacts.items():

            aa1_indices.append(contact['instances'][0]['data']['aa1']['index'])
            plddts = []

            aa1_plddts = []
            aa2_plddts = []
            paes = []
            aa1 = None
            aa2 = None
            models = ''
            for instance in contact['instances']:

                aa1 = instance['data']['aa1'] 
                aa2 = instance['data']['aa2']
                aa1_plddts.append(aa1['plddt'])
                aa2_plddts.append(aa2['plddt'])
                paes.append(instance['pae'])
                plddts.append(mean([aa1['plddt'], aa2['plddt']]))
                models += "," + instance['found_in_model']

            plddt_averages.append(mean(plddts))
        
            line = [
                complex_name,  
                aa1['type'],
                aa1['index'],
                aa1['atom'],
                round(mean(aa1_plddts), 1),
                aa2['type'],
                aa2['index'],
                aa2['atom'],
                round(mean(aa2_plddts), 1),
                contact['instances'][-1]['data']['min_distance'],
                round(mean(paes), 1),
                str(round(100*len(contact['instances'])/len(pdb_file_names), 0)) + '%',
                models
            ]
            output_txt += '\t'.join( [str(i) for i in line]) + '\n'

        num_complete += 1
        print("Finished (" + str(round(100*num_complete/len(complex_names), 1)) + "%) " + complex_name);

        aa1_indices = list(set(aa1_indices))
        aa1_indices.sort()

        interfaces = []
        mean_plddt = mean(plddt_averages)

        if len(aa1_indices) > 0:
            start = min(aa1_indices)
            for i in range(0, len(aa1_indices) -1):
                if(aa1_indices[i + 1] > aa1_indices[i] + 1):

                    if aa1_indices[i] - start > 2:
                        interfaces.append({'start':start, 'end':aa1_indices[i]})
                    
                    start = aa1_indices[i + 1]
            
            if aa1_indices[-1] - start > 2:
                interfaces.append({'start':start, 'end':aa1_indices[-1]})


        pdoq = 0
        if len(aa1_indices) > 1:
            pdoq = pdoq_score(mean_plddt*math.log10(len(aa1_indices)))
        summary_output_data.append({'name': complex_name,'num_aa_contacts':len(aa1_indices), 'interfaces':interfaces, 'mean_plddt':mean_plddt, 'pdoq':pdoq})

        if put_in_folders:
            dest_folder = join(path, complex_name) + '/'
            if os.path.exists(dest_folder) == False:
                 mkdir(dest_folder)
            for filename in get_all_files_for_complex(path, complex_name):
                
                if filename == complex_name:
                    continue
                shutil.move(join(path,filename), dest_folder + filename)
        

    output_file = open(output_filename + ".tsv", "w")
    output_file.write(output_txt)
    output_file.close()

    if len(summary_output_data) > 1:
        summary_output = open(output_filename + "_summary.tsv", "w")
        summary_output.write('name\tnum_aas_in_contacts\tmean_plddt\tpdockq\tmulti_aa_interfaces\n')
        for d in summary_output_data:

            interface_str = ''
            for interface in d['interfaces']:
                interface_str += str(interface['start']) + ':' + str(interface['end']) + ','

            summary_output.write(d['name'] + '\t' + str(d['num_aa_contacts']) + '\t' + str(round(d['mean_plddt'], 1)) + '\t' + str(round(d['pdoq'], 3)) + '\t' + interface_str[:-1] + '\n')

        summary_output.close()
    
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

parser.add_argument(
    "--organize-into-folders",
    action='store_true',
    help="Place all results for a complex into a single folder for that complex",
)

args = parser.parse_args()

process_folder(args.input, args.distance_cutoff, args.pae_cutoff,args.plddt_cutoff,args.output_name, args.organize_into_folders)