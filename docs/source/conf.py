# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(1, os.path.abspath('../../qod_fvm/'))
sys.path.insert(1, os.path.abspath('../../qod_fvm/cell.py'))
sys.path.insert(1, os.path.abspath('../../qod_fvm/field.py'))
sys.path.insert(1, os.path.abspath('../../qod_fvm/mesh.py'))
sys.path.insert(1, os.path.abspath('../../qod_fvm/solver.py'))
sys.path.insert(1, os.path.abspath('../../qod_fvm/data/'))
print(sys.path)


# -- Project information -----------------------------------------------------

project = 'qod_fvm'
copyright = '2020, Quentin Douasbin'
author = 'Quentin Douasbin'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.imgmath', 'recommonmark']
#extensions = ['sphinx.ext.autodoc']
imgmath_image_format = 'svg'
imgmath_font_size = 16

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add suffix to handle markdown
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}
