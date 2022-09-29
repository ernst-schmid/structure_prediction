from argparse import ArgumentParser
import os
import shutil

def get_all_files_for_complex(path, complex_name, filetype = ''):

    files_found = [f for f in os.listdir(path) if complex_name in f and f.endswith(filetype)]
    return files_found

def process_folder(folder_path):

    onlyfiles = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    complex_names = []
    for name in onlyfiles:

        if '.done.txt' in name:
            cname = name[:-9].rstrip("_ -")
            complex_names.append(cname)
    
    output_dir = folder_path + '_analysis'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for complex_name in complex_names:

        complex_directory = os.path.join(output_dir, complex_name)
        if not os.path.exists(complex_directory):
            os.mkdir(complex_directory)

        pdb_filenames = get_all_files_for_complex(folder_path, complex_name, '.pdb')
        for pdb_filename in pdb_filenames:

            new_file_name = pdb_filename.replace("HUMAN", "H").replace("unrelaxed", "ur").replace("rank_", "r").replace("model_", "m")
            shutil.copy(os.path.join(folder_path, pdb_filename), os.path.join(complex_directory, new_file_name))

        img_filenames = get_all_files_for_complex(folder_path, complex_name, '.png')
        for img_filename in img_filenames: 
            new_img_filename = img_filename.replace("HUMAN", "H")
            shutil.copy(os.path.join(folder_path, img_filename), os.path.join(complex_directory, new_img_filename))


parser = ArgumentParser()
parser.add_argument(
    "input",
    default="output_structs",
    help="Directory with structure outputs, pdb_files, from Colabfold",
)

args = parser.parse_args()
process_folder(args.input)