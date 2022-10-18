#!/usr/bin/python
import sys, re, math, numpy as np
from argparse import ArgumentParser
import time 

def dist2(v1, v2):
    return (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2

def atom_from_line(atom_line):

    index = atom_line[23:27]
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


def get_contacts(pdb_filename, a_distance):

    a_distance2 = a_distance**2;    
    threshold = (a_distance + 25)**2
    
    contacts = []
    chains = {}
    residues = {};
    last_chain = None

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
                    chains[chain] = [atom]
                    last_chain = chain

                chains[chain].append(atom)
            
            residues[chain + atom['aa_index']]['atoms'].append(atom)

    num_chains = len(chains)    
    chainnames = list(chains.keys())
    
    for i in range(0, num_chains):
    
        chain1name = chainnames[i]
        chain1 = chains[chain1name]

        for i2 in range(i + 1, num_chains):

            chain2name = chainnames[i2]
            chain2 = chains[chain2name]

            for n1 in chain1:

                rcode1 = chain1name + n1['aa_index']
                r1 = residues[rcode1]
                
                for n2 in chain2:
                    d2 = dist2(n1['xyz'], n2['xyz'])
                    if(d2 < threshold):

                        rcode2 = chain2name + n2['aa_index']
                        r2 = residues[rcode2]
                        min_d2, atoms = get_closest_atoms(r1, r2)
                        if(min_d2 < a_distance2):
                            contacts.append({
                                'distance': round(math.sqrt(min_d2), 1),
                                "aa1":{"chain":chain1name, "type":r1["type"], "index":int(n1["aa_index"]), "abs_index":r1['index'], "atom":atoms[0]['type']},
                                "aa2":{"chain":chain2name, "type":r2["type"], "index":int(n2["aa_index"]), "abs_index":r2['index'], "atom":atoms[1]['type']}
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
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
# print(get_contacts(args.input))


