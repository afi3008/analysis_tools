import re
import os
import glob
import argparse


parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--path_structures",
    type=str,
    required=True,
    help="path to the folder containing the pdb files.",
)
parser_arg.add_argument(
    "--output_path",
    type=str,
    required=True,
    help="path to the pdb file containing the combined pdb",
)


def return_pattern(filepath):
    """
    Returns the index number of the file name.
    :param filepath: str, path of the file
    :return: integer, the index of the file.
    """
    if not isinstance(filepath, str):
        raise TypeError("The file path must of type str")
    if not re.match(".*structure_z_[0-9]+\\.pdb", filepath):
        raise ValueError(
            "The file names must match the format '*/structure_z_[0-9]+\\.pdb'"
        )
    file_name = re.search(r"structure_z_\d+\.pdb", filepath).group(0)
    return int(re.search("[0-9]+", file_name).group(0))


def get_indexes(list_file):
    """
    Get the integer part of the filenames, where the filenames are in the format:
    structure_z_i.pdb where i in the integer.
    :param list_file: list of string, path of each file.
    :return: list of integer, containing the integer of each filename in the same order.
    """
    index_list = [return_pattern(f) for f in list_file]
    return index_list


def read_file(file_path):
    """
    Read the content of a PDB file
    :param file_path: str, path to the file to read.
    :return: str, content of the pdb file.
    """
    with open(file_path, "r") as f:
        content = f.readlines()

    return content


def write_file(content, file_path, i, is_last):
    """
    Write the content in the given file.
    :param content: list of string, lines to write in the file.
    :param file_path: str, path to the file in which we want to write.
    :param i: number of the model.
    :return: None
    """
    with open(file_path, "a") as f:
        f.write(f"MODEL      {i}\n")
        f.writelines(content)
        f.write("ENDMDL\n")
        if is_last:
            f.write("END\n")


def get_files(path_to_structures):
    """
    Get the list of the files in the structure folder and return it
    :param path_to_structures: str, path to the folder containing the structures.
    :return: list of the files path
    """
    if not isinstance(path_to_structures, str):
        raise ValueError(
            "The path to the folder containing the structures should be a string."
        )
    if not os.path.isdir(path_to_structures):
        raise ValueError(f"This path {path_to_structures} does not exists")

    path_format = f"{path_to_structures}structure_z_[0-9]*.pdb"
    structures_path = glob.glob(path_format)
    if structures_path == []:
        raise Exception(f"""The folder {path_to_structures} does not contain any file in the format
                        {path_format}""")
    return structures_path


def sort_files(list_files, list_indexes):
    """
    Sort the list of files by the number of the structure. The file names are in the format:
    structure_z_i.pdb where i is an integer.
    :param list_files: list of str, list of path to the pdb.
    :param list_indexes: list of integer, list of the pdb numbers.
    :return: list of str, list of path to the pdb but sorted by their structure number.
    """
    if not isinstance(list_files, list):
        raise TypeError("Firt argument must be a list")

    if not isinstance(list_indexes, list):
        raise TypeError("Second argument must be a list")

    index_file_tuple = zip(list_indexes, list_files)
    sorted_index_file_tuple = sorted(index_file_tuple)
    sorted_indexes, sorted_files = list(zip(*sorted_index_file_tuple))
    return sorted_indexes, sorted_files


def check_output_file_exists(path):
    """
    Check if the file exists and deletes it if it does
    :param path: str, path to the file that will combine the different pdb.
    :return: None
    """
    if os.path.isfile(path):
        os.remove(path)


def combine(path_to_structures, path_to_combined_file):
    """
    This function combines different pdb files in a single one
    containing different models.
    Arguments:
        - path_to_structures: str, path to the folder containing the structures
        - path_to_combined_file: str, path to the file that will combine the contents.
    :return: None
    """
    _, ext = os.path.splitext(path_to_combined_file)
    if ext != ".pdb":
        raise Exception("The output file must be a .pdb file.")

    list_files = get_files(path_to_structures)
    list_indexes = get_indexes(list_files)
    sorted_list_indexes, sorted_list_files = sort_files(list_files, list_indexes)
    N_files = len(sorted_list_files)
    check_output_file_exists(path_to_combined_file)
    for i, (index, file) in enumerate(zip(sorted_list_indexes, sorted_list_files)):
        is_last = N_files == i + 1
        content = read_file(file)
        write_file(content, path_to_combined_file, index, is_last)


if __name__ == "__main__":
    args = parser_arg.parse_args()
    struct_path = args.path_structures
    output_path = args.output_path
    combine(struct_path, output_path)
