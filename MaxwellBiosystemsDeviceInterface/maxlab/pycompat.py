import sys

def is_py2():
    return True if sys.version_info[0] == 2 else False

def decode(obj):
    return obj if is_py2() else obj.decode('latin1')

def encode(obj):
    return obj if is_py2() else str.encode(obj, 'latin1')


