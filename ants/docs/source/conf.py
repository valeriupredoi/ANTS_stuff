# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
#
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from datetime import datetime

sys.path.insert(0, os.path.abspath("../../bin"))  # For autodoc
sys.path.insert(0, os.path.abspath("../../lib"))  # For autodoc

import ants  # noqa: E402

# -- Project information -----------------------------------------------------

project = "ANTS"
copyright = "(C) British Crown Copyright 2015 - {}, Met Office".format(
    datetime.now().year
)
author = "Harold Dyson, Esther Turner, Katherine Tomkins, Andrew Clark"

# The full version, including alpha/beta/rc tags
version = ants.__version__
release = ants.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    "sphinx.ext.viewcode",
    "sphinxarg.ext",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

intersphinx_mapping = {
    "python": (
        "https://docs.python.org/{}.{}".format(
            sys.version_info.major, sys.version_info.minor
        ),
        None,
    ),
    "iris": ("https://scitools-iris.readthedocs.io/en/v3.0.0", None),
    "mule": ("https://code.metoffice.gov.uk/doc/um/mule/latest", "mule-objects.inv"),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

extlinks = {
    "anciltrac": ("https://code.metoffice.gov.uk/trac/ancil/%s", "ancil trac %s"),
    "ancilwiki": ("https://code.metoffice.gov.uk/trac/ancil/wiki/%s", "ancil wiki %s"),
    "contrib": (
        "https://code.metoffice.gov.uk/trac/ancil/browser/contrib/trunk/%s",
        "contrib %s",
    ),
    "fcm": ("http://metomi.github.io/fcm/doc/%s", "fcm %s"),
    "source": (
        "https://code.metoffice.gov.uk/trac/ancil/browser/ants/trunk/%s",
        "source %s",
    ),
}

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "special-members": "__init__, __call__",
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for link checking -----------------------------------------------

linkcheck_ignore = [
    r"https://code.metoffice.gov.uk/trac/ancil/newticket.*",
    # Needed due to authentication and redirects on MOSRS breaking linkcheck
    # when linking to specific anchors.  The URLs that need to be matched are
    # currently:
    # https://code.metoffice.gov.uk/doc/um/mule/latest/mule.html#mule.UMFile
    # https://code.metoffice.gov.uk/doc/um/mule/latest/mule.html#mule.Field3
    # https://code.metoffice.gov.uk/doc/um/mule/latest/mule/ancil.html#mule.ancil.AncilFile
    r"https://code.metoffice.gov.uk/doc/um/mule/latest/mule.*\.html\#",
]
