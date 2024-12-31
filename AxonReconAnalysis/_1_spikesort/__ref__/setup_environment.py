import os
import sys

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

def set_pythonpath():
    #workspace_folder = os.path.dirname(os.path.abspath(__file__))
    workspace_folder = get_git_root_ws()
    os.environ['PYTHONPATH'] = workspace_folder
    sys.path.insert(0, workspace_folder)

if __name__ == "__main__":
    set_pythonpath()