import requests, pandas as pd
import urllib

from argparse import ArgumentParser
from pdb_parser import atom_from_line

def get_plddts(entry):

    url_start = 'https://alphafold.ebi.ac.uk/files/AF-'
    url_end = '-F1-model_v3.pdb'
    pdb_url = url_start + entry + url_end
    plddts = []

    try:
        pdb_file = urllib.request.urlopen(pdb_url)
        for line in pdb_file:

            atom_line = line.decode("utf-8")

            if('ATOM' not in atom_line):
                continue

            atom = atom_from_line(atom_line)
            if(atom['type'] == 'N'):
                plddts.append(atom['bfactor'])
    except urllib.error.HTTPError as err:
        print(err)

    return plddts
        

parser = ArgumentParser()
parser.add_argument(
    "input",
    default="entries.tsv",
    help="TSV file with list of Uniprot entries",
)
args = parser.parse_args()

entries_df = pd.read_csv(args.input, delimiter='\t')

output_txt = 'entry_name\tentry\t'
output_filename = 'alphafold_plddts'
max_num_aas = entries_df['Length'].max()
for i in range(1, max_num_aas + 1):
    output_txt += str(i) + '\t'

output_txt += '\n'

for index, row in entries_df.iterrows():

    entry = row['Entry']
    entry_name = row['Entry Name']
    print('started ' + entry_name)
    plddts = get_plddts(row['Entry'])
    output_txt += entry_name + '\t' + entry + '\t' + '\t'.join([str(p) for p in plddts]) + '\n'
    print('finished ' + entry_name)

output_file = open(output_filename + ".tsv", "w")
output_file.write(output_txt)
output_file.close()
