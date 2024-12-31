import os
import sys
from pprint import pprint

#use git to get root directory of repository
def get_git_root():
    git_root = os.popen('git rev-parse --show-toplevel').read().strip()
    return git_root

def get_git_root_ws():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    while current_dir != os.path.dirname(current_dir):  # Root directory check
        if os.path.isdir(os.path.join(current_dir, '.git')):
            return current_dir
        current_dir = os.path.dirname(current_dir)
    return None

def add_external_paths():
    workspace_folder = get_git_root_ws()
    external_folder = os.path.join(workspace_folder, 'external')
    for sub_external_folder in os.scandir(external_folder):
        if sub_external_folder.is_dir():
            sys.path.insert(0, sub_external_folder.path)

def set_pythonpath():
    #pprint(sys.path)
    workspace_folder = get_git_root_ws()
    sys.path.insert(0, workspace_folder)
    add_external_paths()
    #pprint(sys.path)

if __name__ == "__main__":
    set_pythonpath()