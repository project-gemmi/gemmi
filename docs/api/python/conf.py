PROJECT_TITLE = "Gemmi"
PROJECT_SUBTITLE = "Python API"

INPUT_MODULES = ['gemmi']

MAIN_PROJECT_URL = "../.."
PYBIND11_COMPATIBILITY = True


PLUGINS = ['m.sphinx']

M_SPHINX_INVENTORIES = [...]
M_SPHINX_INVENTORY_OUTPUT = 'gemmi.inv'
M_SPHINX_PARSE_DOCSTRINGS = False

M_SPHINX_INVENTORIES = [
    ('sphinx/python.inv', 'https://docs.python.org/3/', ['xml.']),
    ('sphinx/numpy.inv', 'https://docs.scipy.org/doc/numpy/', [''], ['m-flat'])
]

LINKS_NAVBAR1 = [
    ('Modules', 'modules', []),
    ('Classes', 'classes', []),
    ('Github', 'https://github.com/project-gemmi/gemmi', [])
]

OUTPUT = "../../_build/html/api/python"