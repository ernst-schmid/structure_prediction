#!/usr/bin/python
import sys, re, math
from argparse import ArgumentParser
import time 

def dist2(v1, v2):
    return (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2

def atom_from_line(atom_line):

    chain = atom_line[20:22].strip()
    index = int(atom_line[23:27])
    coordinates = [float(i) for i in re.findall(r"[-+]?(?:\d*\.\d+|\d+)", atom_line[29:54])]
    aa = atom_line[17:20]
    bfactor = float(atom_line[60:66])

    return {"type":atom_line[13:16].strip(), "chain":chain, "aa_index":index, "aa":aa, "coordinates":coordinates, "bfactor":bfactor}

def get_closest_atoms(res1, res2):

    min_d2 = 1e6
    atoms = [None, None]
    for a1 in res1:
         for a2 in res2:
            d2 = dist2(a1["coordinates"], a2["coordinates"])
            if d2 < min_d2:
                min_d2 = d2
                atoms[0] = a1
                atoms[1] = a2
    return atoms


def get_contacts(pdb_filename, a_distance):

    a_distance2 = a_distance**2;    
    threshold = (a_distance + 30)**2
    
    contacts = []
    chains = {}
    residues = {};

    with open(pdb_filename) as pdb_file:

        res_num = 0
        for atom_line in pdb_file:

            if('ATOM' not in atom_line):
                continue
            
            atom = atom_from_line(atom_line)
            if(atom['type'] == 'N'):
                res_num += 1
            
            atom['absolute_aa_index'] = res_num

            rescode = atom['chain'] + ':' + atom['aa'] + str(atom['aa_index'])
            if(rescode in residues):
                residues[rescode].append(atom)
            else:
                residues[rescode] = [atom]

            if('CA' in atom['type']):
                if atom['chain'] in chains:
                    chains[atom['chain']].append(atom)
                else:
                    chains[atom['chain']] = [atom]

    num_chains = len(chains)    
    chainnames = list(chains.keys())
    
    for i in range(0, num_chains):
    
        chain1 = chains[chainnames[i]]

        for i2 in range(i + 1, num_chains):
            chain2 = chains[chainnames[i2]]

            for ca1 in chain1:
                
                for ca2 in chain2:
                    d2 = dist2(ca1['coordinates'], ca2['coordinates'])
                    if(d2 < threshold):

                        rescode1 = ca1['chain'] + ':' + ca1['aa'] + str(ca1['aa_index'])
                        rescode2 = ca2['chain'] + ':' + ca2['aa'] + str(ca2['aa_index'])
                        atoms = get_closest_atoms(residues[rescode1], residues[rescode2])

                        min_distance = dist2(atoms[0]['coordinates'], atoms[1]['coordinates'])
                        if(min_distance < a_distance2):
                            contacts.append({
                                "ca_distance": round(math.sqrt(d2), 1),
                                'min_distance': round(math.sqrt(min_distance), 1),
                                "aa1":{"chain":chainnames[i], "type":ca1["aa"], "index":ca1["aa_index"],  "abs_index":ca1["absolute_aa_index"], "atom":atoms[0]['type']},
                                "aa2":{"chain":chainnames[i2], "type":ca2["aa"], "index":ca2["aa_index"], "abs_index":ca2["absolute_aa_index"], "atom":atoms[1]['type']}
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