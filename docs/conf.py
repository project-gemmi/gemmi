# -*- coding: utf-8 -*-

# -- General configuration ------------------------------------------------

# while we use Sphinx 8+, old version suffices to run doctests
needs_sphinx = '5.3.0'

extensions = ['sphinx.ext.doctest']

templates_path = ['_templates']

master_doc = 'index'

project = u'Gemmi'
copyright = u'Global Phasing Ltd'
author = u'Marcin Wojdyr'

with open('../include/gemmi/version.hpp') as _f:
    for _line in _f:
        if _line.startswith('#define GEMMI_VERSION '):
            version = _line.split()[2].strip('"')
release = version

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
highlight_language = 'cpp'
default_role = 'literal'


# -- Options for HTML output ----------------------------------------------

html_theme = 'furo'
html_theme_options = {
    "source_repository": "https://github.com/project-gemmi/gemmi/",
    "source_branch": "master",
    "source_directory": "docs/",
}
html_static_path = ['_static']
html_css_files = ['custom.css']

# Edit link can be also used to see the source
html_show_sourcelink = False
html_copy_source = False

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # 'papersize': 'letterpaper',
    # 'pointsize': '10pt',
    # 'preamble': '',
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'Gemmi.tex', u'Gemmi Documentation',
     u'Marcin Wojdyr', 'manual'),
]

doctest_global_setup = '''
import os
import sys
disabled_features = []
try:
    import numpy
    if numpy.__version__ >= '2':
        numpy.set_printoptions(legacy='1.25')
    numpy.set_printoptions(threshold=5)
except ImportError:
    disabled_features.append('NumPy')
    numpy = None
try:
    import pandas
except ImportError:
    disabled_features.append('pandas')
    pandas = None
try:
    import networkx
except ImportError:
    disabled_features.append('networkx')
    networkx = None
try:
    import pynauty
except ImportError:
    disabled_features.append('pynauty')
    pynauty = None
ccp4_path = os.getenv('CCP4')
if ccp4_path:
    mdm2_unmerged_mtz_path = (ccp4_path + '/lib/python3.9/site-packages/'
                              + 'ccp4i2/demo_data/mdm2/mdm2_unmerged.mtz')
    if not os.path.isfile(mdm2_unmerged_mtz_path):
        ccp4_path = None
if ccp4_path is None:
    disabled_features.append('$CCP4')

import gemmi
gemmi.set_leak_warnings(False)
'''
