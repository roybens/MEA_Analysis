import json
import numpy as np
import math
import os
import fnmatch

def get_surrounding_coordinates(x, y, number):
    surrounding_coordinates = []
    number = math.sqrt(number)
    for i in range(-2, 3):
        for j in range(-2, 3):
            surrounding_coordinates.append((x + i, y + j))
    
    return surrounding_coordinates


def convert_int64_keys_to_ints(d):
    """
    Recursively convert all numpy.int64 keys in a dictionary to integers.
    """
    new_dict = {}
    for k, v in d.items():
        if isinstance(k, np.int64):
            k = int(k)
        if isinstance(v, dict):
            v = convert_int64_keys_to_ints(v)
        new_dict[k] = v
    return new_dict

def load_json(filename):
    d = {}
    with open(filename,'rb') as f:
        d = json.load(f)
    return d

def dumpdicttofile(data,filename):
    data = convert_int64_keys_to_ints(data)
    json_data = json.dumps(data,indent=4)

    with open(filename,'w') as fileptr:
        fileptr.write(json_data)

def get_key_by_value(d, value):
    """
    Useful to get the templates associated with a extreme electrode
    """
    keys = [k for k, v in d.items() if v == value]
    return keys if keys else None



def inarrxnotinarry(arr1,arr2):
    print(f" in array1 not in array2 :{[x for x in arr1 if x not in arr2]} ")
    print(f" in array2 not in array1 :{[x for x in arr2 if x not in arr1]} ")


def flatten(lst):
    flat_list = []
    for item in lst:
        if isinstance(item, list):
            flat_list.extend(flatten(item))
        else:
            flat_list.append(item)
    return flat_list

def get_templates_with_same_channels(electrode_file):

    my_dict = load_json(electrode_file)

    templates_with_channel = {}

    for template_name, template_data in my_dict.items():
        for channel_name in template_data:
            if channel_name in templates_with_channel:
                templates_with_channel[channel_name].append(template_name)
            else:
                templates_with_channel[channel_name] = [template_name]
    
    same_channel_templates = [templates for channel, templates in templates_with_channel.items() if len(templates) > 1]

    return same_channel_templates


def empty_directory(directory_path):
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        # List all files and subdirectories in the directory
        file_list = os.listdir(directory_path)
        
        for filename in file_list:
            file_path = os.path.join(directory_path, filename)
            if os.path.isfile(file_path):
                # Remove files
                os.remove(file_path)
            elif os.path.isdir(file_path):
                # Remove subdirectories and their contents recursively
                empty_directory(file_path)
                os.rmdir(file_path)
        print(f"The directory '{directory_path}' has been emptied.")
    else:
        print(f"The directory '{directory_path}' does not exist.")



def find_files_with_subfolder(root_dir, file_name_pattern, subfolder_name):
    file_paths = []
    for dirpath, _, filenames in os.walk(root_dir):
        if subfolder_name in dirpath.split(os.path.sep):
            for filename in fnmatch.filter(filenames, file_name_pattern):
                file_paths.append(os.path.join(dirpath, filename))
    return file_paths

def isexists_folder_not_empty(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print("sorting dir doesnt exist: making new one")
        os.makedirs(folder_path,0o777)
        return 0
    # Get the list of items (files and folders) in the specified directory
    items = os.listdir(folder_path)

    # If the folder is not empty, return True
    return len(items) > 0


