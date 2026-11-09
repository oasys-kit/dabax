# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import dabax

# -- Project information -----------------------------------------------------

project = 'dabax'
copyright = '2024-2026, M. Sanchez del Rio'
author = 'M. Sanchez del Rio'
import tomllib
with open(os.path.join(os.path.dirname(__file__), '..', 'pyproject.toml'), 'rb') as _f:
    release = tomllib.load(_f)['project']['version']
version = release

# -- General configuration ---------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx.ext.autosummary',
              'myst_parser'
]

myst_enable_extensions = ["dollarmath"]

autosummary_generate = True

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_context = {
    'display_github': True,
    'github_user': 'oasys-kit',
    'github_repo': 'dabax',
    'github_version': 'main/docs/',
}


# The pages in generated/ are created at build time by autosummary and the index pages
# have no source file, so "Edit on GitHub" would give a 404 for them.
def _hide_github_link(app, pagename, templatename, context, doctree):
    if pagename.startswith('generated/') or pagename in ('genindex', 'py-modindex', 'search'):
        context['display_github'] = False


def setup(app):
    app.connect('html-page-context', _hide_github_link)
