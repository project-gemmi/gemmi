import os
import tempfile

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

def get_path_for_tempfile(suffix=''):
    handle, out_name = tempfile.mkstemp(suffix=suffix)
    os.close(handle)
    return out_name
