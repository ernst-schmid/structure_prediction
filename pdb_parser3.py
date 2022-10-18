#!/usr/bin/python
import sys, re, math, numpy as np
from argparse import ArgumentParser
import time 

def dist2(v1, v2):
    return (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2

def atom_from_line(atom_line):

    index = atom_line[23:27].strip()
    coordinates = [float(i) for i in re.findall(r"[-+]?(?:\d*\.\d+|\d+)", atom_line[29:54])]

    return {"type":atom_line[13:16].strip(), "aa_index":index, "xyz":coordinates,}

def get_closest_atoms(res1, res2):

    min_d2 = 1e6
    atoms = [None, None]
    for a1 in res1['atoms']:
         for a2 in res2['atoms']:
            d2 = dist2(a1["xyz"], a2["xyz"])
            if d2 < min_d2:
                min_d2 = d2
                atoms[0] = a1
                atoms[1] = a2
    return (min_d2, atoms)


def get_closest_atoms2(res1, res2):

    r1_coords = []
    for a in res1['atoms']:
        r1_coords.append(a['xyz'])

    r1_coords = np.array(r1_coords)
    r1_coords = np.repeat(r1_coords[None, ...],len(res2['atoms']), 0)

    r2_coords = []
    for a in res2['atoms']:
        r2_coords.append(a['xyz'])

    r2_coords = np.array(r2_coords)
    r2_coords = np.tile(r2_coords, (1, len(res1['atoms']))).reshape(len(res2['atoms']),len(res1['atoms']),3)
    d2s = np.sum(((r1_coords - r2_coords)**2), 2)
    min_d2 = d2s.min()
    y, x = np.where(d2s == min_d2)

    return (min_d2, [res1['atoms'][x[0]], res2['atoms'][y[0]]])


def get_contacts(pdb_filename, a_distance):

    a_distance2 = a_distance**2;    
    threshold = (a_distance + 25)**2
    
    contacts = []
    chainnames = []
    residues = {};
    last_chain = None
    chain_index = -1

    coordinates = []

    with open(pdb_filename) as pdb_file:

        res_num = 0
        for atom_line in pdb_file:

            if (atom_line[0:4] != 'ATOM'):
                continue
            
            atom = atom_from_line(atom_line)
            if(atom['type'] == 'N'):

                chain = atom_line[20:22].strip()
                res_num += 1
                residues[chain + atom['aa_index']] = {"atoms":[], "index":res_num, "type": atom_line[17:20]}

                if chain != last_chain:
                    chainnames.append(chain)
                    last_chain = chain
                    chain_index += 1
                    coordinates.append([])

                coordinates[chain_index].append(atom['xyz'])
            
            residues[chain + atom['aa_index']]['atoms'].append(atom)

    num_chains = chain_index  

    for i in range(0, num_chains):
        coordinates[i] = np.array(coordinates[i])
    
    for i in range(0, num_chains):
    
        chain1name = chainnames[i]
        c1_coords = coordinates[i]

        for i2 in range(i + 1, num_chains):

            chain2name = chainnames[i2]
            c2_coords = coordinates[i2]

            for aa1_index in range(0, len(c1_coords)):

                r1 = residues[chain1name + str(aa1_index + 1)]
                d2s = ((np.tile(c1_coords[aa1_index], (len(c2_coords), 1)) - c2_coords)**2).sum(axis=1)
                aa2_indices = np.where(d2s < threshold)[0]
                for aa2_index in aa2_indices:

                    r2 = residues[chain2name + str(aa2_index + 1)]
                    min_d2, atoms = get_closest_atoms2(r1, r2)
                    if(min_d2 < a_distance2):
                        contacts.append({
                            'distance': round(math.sqrt(min_d2), 1),
                            "aa1":{"chain":chain1name, "type":r1["type"], "index":aa1_index, "abs_index":r1['index'], "atom":atoms[0]['type']},
                            "aa2":{"chain":chain2name, "type":r2["type"], "index":aa2_index, "abs_index":r2['index'], "atom":atoms[1]['type']}
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
get_contacts(args.input, 10)
# print(get_contacts(args.input, 3))
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))

