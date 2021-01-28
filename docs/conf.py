# -*- coding: utf-8 -*-
import os

# -- General configuration ------------------------------------------------

needs_sphinx = '1.8.5'

extensions = ['sphinx.ext.doctest', 'sphinx.ext.githubpages']

#templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

project = u'Gemmi'
copyright = u'2017-2021 Global Phasing Ltd'
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

if not os.environ.get('READTHEDOCS'):
    html_theme = 'sphinx_rtd_theme'
    #import cloud_sptheme as csp
    #html_theme = "cloud"
    #html_theme_path = [csp.get_theme_dir()]
    #html_theme = 'bizstyle'
    # html_theme_options = {}

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
import sys
assert sys.version_info[0] > 2, "Tests in docs are for Python 3 only"
try:
    import numpy
except ImportError:
    print('Tests that use NumPy are disabled.', file=sys.stderr)
    numpy = None
try:
    import pandas
except ImportError:
    print('Tests that use pandas are disabled.', file=sys.stderr)
    pandas = None
try:
    import networkx
except ImportError:
    print('Tests that use networkx are disabled.', file=sys.stderr)
    networkx = None
import os
mdm2_unmerged_mtz_path = os.getenv('CCP4')
if mdm2_unmerged_mtz_path:
    mdm2_unmerged_mtz_path += '/share/ccp4i2/demo_data/mdm2/mdm2_unmerged.mtz'
    if not os.path.isfile(mdm2_unmerged_mtz_path):
        mdm2_unmerged_mtz_path = None
'''

def setup(app):
    app.add_css_file('custom.css')
