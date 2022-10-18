#!/usr/bin/python
import numpy as np
import re
from argparse import ArgumentParser

def atom_from_line(atom_line):

    chain = atom_line[20:22].strip()
    index = int(atom_line[23:27])
    coordinates = np.array([float(i) for i in re.findall(r"[-+]?(?:\d*\.\d+|\d+)", atom_line[29:54])])
    aa = atom_line[17:20]
    return {"type":atom_line[13:16].strip(), "chain":chain, "idx":index, "aa":aa, "xyz":coordinates}

def dist2(xyz1, xyz2):

    return np.sum((xyz1-xyz2)**2, axis=0)

def atom_dist2(a1, a2):

    return np.sum((a1["xyz"]-a2["xyz"])**2, axis=0)

def get_closest_atoms(res1, res2):

    min_d2 = 1e6
    atoms = [None, None]
    for a1 in res1['atoms']:
         for a2 in res2['atoms']:
            d2 = atom_dist2(a1, a2)
            if d2 < min_d2:
                min_d2 = d2
                atoms[0] = a1
                atoms[1] = a2
    return atoms

def get_contacts(pdb_filename, a_dist = 5):

    residues = {}
    last_res_chain = None
    chainnames = [];

    d2_large = (a_dist + 20)**2

    coords = []
    last_chain = None
    chain_index = -1
    with open(pdb_filename) as pdb_file:

        for atom_line in pdb_file:
            if (atom_line[0:4] != 'ATOM' or atom_line[13:16].strip() != 'N'):
                continue
            
            chain = atom_line[20:22].strip()
            if chain != last_chain:
                last_chain = chain
                chain_index += 1
                coords.append([])

            xyz = np.array([float(i) for i in re.findall(r"[-+]?(?:\d*\.\d+|\d+)", atom_line[29:54])])
            coords[chain_index].append(xyz)
    

    rough_contacts = []

    for c1 in range(0, chain_index + 1):

        c1_coords = coords[c1]
        for c2 in range(c1 + 1, chain_index + 1):
            
            c2_coords = coords[c2]
            for a1_idx, a1 in enumerate(c1_coords):
                 for a2_idx, a2 in enumerate(c2_coords):
                    if dist2(a1, a2) < d2_large:
                        rough_contacts.append([c1, a1_idx, c2, a2_idx])
    
    return rough_contacts


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="struct_path",
    help="PDB file of structure",
)
args = parser.parse_args()
get_contacts(args.input)
# print(get_contacts(args.input))



