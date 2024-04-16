"""
quick_install_editable.py

This script is used to quickly update and install the requirements associated with the MEAProcessingLibrary 
during development.

The script uses pipreqs, a package that automatically updates requirements.txt based on the imports in your 
project. It is required to be installed in order to run this script.

The script also uses 'pip install -e' to install MEAProcessingLibrary in editable mode. 
This means that changes to the code will be reflected in the package without needing to reinstall it. 
This works by creating a symbolic link to the package in the site-packages folder of your Python environment.
"""

# Standard library imports
import subprocess
import os

try:   
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Run pipreqs and overwrite requirements.txt
    subprocess.run(['pipreqs', script_dir, '--force'])

    # Install the current package in editable mode
    subprocess.run(['pip', 'install', '-e', script_dir])
except Exception as e:
    print("Failed to update requirements.txt and install package.")
    print("Error: ", e)
    print("Make sure you have the current version of pip installed and in your path.")
    print("Then, manually run the following command in terminal: 'pip install pipreqs'")
    print("Then, try running this quick install again.")