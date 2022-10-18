from hashlib import new
import os
import pandas as pd
from argparse import ArgumentParser

def process_file(rename_file):
    
    csv_data = pd.read_csv(rename_file)
    for index, row in csv_data.iterrows():

        if index > 491:
            break

        old_str = row['old']
        h = round(0.5*len(row['new']))
        mix_str = '$$$$'
        new_str = row['new'][0:h] + mix_str + row['new'][h:]
        names = [];
        sub_df = csv_data.iloc[index:, :]
        for index, row in sub_df.iterrows():
            if(old_str in row['old']):
                names.append(row['old'].replace(old_str, new_str))
        
        if len(names) > 0:
            print(name)
            print(names)
        

parser = ArgumentParser()
parser.add_argument(
    "--rename-file",
    help="CSV file of renaming map",
    type=str,
)
args = parser.parse_args()
process_file( args.rename_file)