from hashlib import new
import os
import pandas as pd
from argparse import ArgumentParser
import re

def process_file(folder):
    
    all_files = os.listdir(folder)
    
    for f in all_files:
        new_name = re.sub(r"_\d+aa$", "", f) 
        os.rename(f, new_name)


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="predictions",
    help="folder containing files to rename",
)

args = parser.parse_args()
process_file(args.input)