# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
CLASS Reader package, migrated from pyspeckit
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from core import ClassObject
