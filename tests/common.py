import os
import tempfile
try:
    import numpy
    numpy_version = tuple(int(n) for n in numpy.__version__.split('.')[:2])
except ImportError:
    numpy = None

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

def get_path_for_tempfile(suffix=''):
    handle, out_name = tempfile.mkstemp(suffix=suffix)
    os.close(handle)
    return out_name

def assert_numpy_equal(self, arr1, arr2):
    # equal_nan arg was added in NumPy 1.19.0
    if numpy and numpy_version >= (1,19):
        self.assertTrue(numpy.array_equal(arr1, arr2, equal_nan=True))
