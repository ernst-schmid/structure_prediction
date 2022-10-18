#!/usr/bin/python
import numpy as np
import scipy as sp
import re
from argparse import ArgumentParser
import time 

def atom_from_line(atom_line):

    chain = atom_line[20:22].strip()
    index = int(atom_line[23:27])
    coordinates = np.array([float(i) for i in re.findall(r"[-+]?(?:\d*\.\d+|\d+)", atom_line[29:54])])
    aa = atom_line[17:20]
    # bfactor = float(atom_line[60:66])

    return {"type":atom_line[13:16].strip(), "chain":chain, "idx":index, "aa":aa, "xyz":coordinates}

def atom_dist(a1, a2):

    return np.linalg.norm(a1["xyz"] - a2["xyz"])

def get_closest_atoms(res1, res2):

    min_d2 = 1e6
    atoms = [None, None]
    for a1 in res1['atoms']:
         for a2 in res2['atoms']:
            d2 = atom_dist(a1, a2)
            if d2 < min_d2:
                min_d2 = d2
                atoms[0] = a1
                atoms[1] = a2
    return atoms


def get_contacts(pdb_filename, a_dist = 5):
    
    residues = {}
    coords = []
    chainnames = []
    last_ca_chain = None
    chain_index = -1;
    gross_dist = a_dist + 30

    abs_res_idx = 0
    with open(pdb_filename) as pdb_file:

        for atom_line in pdb_file:

            if('ATOM' not in atom_line):
                continue

            atom = atom_from_line(atom_line)
            chain = atom['chain']
            atype = atom['type']
            rescode = chain + ':' + str(atom['idx'])
            if(atype == 'N'):
                residues[rescode] = {"atoms":[atom], "index":abs_res_idx}
                abs_res_idx += 1
            else:
                residues[rescode]['atoms'].append(atom)

            if(atype == 'CA'):
                if chain != last_ca_chain:
                    chain_index += 1
                    coords.append([])
                    chainnames.append(chain)

                coords[chain_index].append(atom['xyz'])
                last_ca_chain = chain
    
    contacts = []
    num_chains = len(coords)    
    for i in range(0, num_chains):

        chain1_name = chainnames[i]

        for i2 in range(i + 1, num_chains):

            chain2_name = chainnames[i2]

            distances = sp.spatial.distance.cdist(coords[i], coords[i2])
            aa1_ixs, aa2_ixs = np.where(distances < gross_dist)
            for cindex in range(0, len(aa1_ixs)):
                res1 = residues[chain1_name + ':' + str(aa1_ixs[cindex] + 1)]
                res2 = residues[chain2_name + ':' + str(aa2_ixs[cindex] + 1)]
                atoms = get_closest_atoms(res1, res2)
                ad = atom_dist(atoms[0], atoms[1])
                if ad <= a_dist:
                    contacts.append({
                            'distance': round(ad, 1),
                            "aa1":{"chain":chain1_name, "type":atoms[0]["aa"], "index":atoms[0]["idx"], "abs_index":res1['index'], "atom":atoms[0]['type']},
                            "aa2":{"chain":chain2_name, "type":atoms[1]["aa"], "index":atoms[1]["idx"], "abs_index":res2['index'], "atom":atoms[1]['type']}
                        })

    return contacts


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="struct_path",
    help="PDB file of structure",
)
args = parser.parse_args()
startTime = time.time()
get_contacts(args.input)
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
