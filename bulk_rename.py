from hashlib import new
import os
import pandas as pd
from argparse import ArgumentParser

def process_file(folder, rename_file):
    
    csv_data = pd.read_csv(rename_file)
    all_files = os.listdir(folder)

    mix_str = '$$$$'
    for index, row in csv_data.iterrows():
        
        old_str = row['old']
        new_str = row['new']
        h = round(0.5*len(row['new']))
        new_str = row['new'][0:h] + mix_str + row['new'][h:]
        for f in all_files:
            if old_str in f:
                new_filename = f.replace(old_str, new_str)
                os.rename(f, new_filename)
    
    for f in all_files:
        os.rename(f, f.replace(mix_str, ''))


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="predictions",
    help="folder containing files to rename",
)
parser.add_argument(
    "--rename-file",
    help="CSV file of renaming map",
    type=str,
)


args = parser.parse_args()
process_file(args.input, args.rename_file)