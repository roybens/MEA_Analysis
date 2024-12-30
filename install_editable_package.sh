#get location of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#set up directory for the package
echo "directory $DIR"

#install the package in editable mode
pip install -e $DIR