# -*- coding: utf-8 -*-

# -- General configuration ------------------------------------------------

needs_sphinx = '5.3.0'

extensions = ['sphinx.ext.doctest', 'sphinx.ext.githubpages']

#templates_path = ['_templates']

source_suffix = '.rst'

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
highlight_language = 'c++'


# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['custom.css']

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
assert sys.version_info[0] > 2, "Tests in docs are for Python 3 only"
disabled_features = []
try:
    import numpy
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
mdm2_unmerged_mtz_path = os.getenv('CCP4')
if mdm2_unmerged_mtz_path:
    mdm2_unmerged_mtz_path += ('/lib/python3.7/site-packages/' +
                               'ccp4i2/demo_data/mdm2/mdm2_unmerged.mtz')
    if not os.path.isfile(mdm2_unmerged_mtz_path):
        mdm2_unmerged_mtz_path = None
if mdm2_unmerged_mtz_path is None:
    disabled_features.append('$CCP4')
'''

def setup(app):
    app.add_css_file('custom.css')
