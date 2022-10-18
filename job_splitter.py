import pandas as pd
from argparse import ArgumentParser

def process_file(file, num_split, max_aa_length):
    
    csv_data = pd.read_csv(file)
    csv_data.sort_values(by=['length'])
    csv_data.drop(csv_data[csv_data['length'] > max_aa_length].index, inplace = True);

    for i in range(1, num_split + 1):

        extracted_file = csv_data.iloc[i::num_split, :]
        extracted_file.to_csv(file + '.split_' + str(i) + '.csv')


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="output_structs",
    help="CSV file of sequences to be split",
)
parser.add_argument(
    "--num-split",
    default=4,
    help="Number of files to generate",
    type=float,
)

parser.add_argument(
    "--max-length",
    default=2350,
    help="Maximum length of complex to be placed into new CSV file",
    type=float,
)


args = parser.parse_args()
process_file(args.input, args.num_split, args.max_length)