"""
adapted from numpy.oldnumeric version 1.9.0 (dropped from numpy after that). 
"""
from __future__ import division, absolute_import, print_function

from .functions import *
from .precision import *
from .ufuncs import *

from . import functions
from . import precision
from . import ufuncs

__all__ = []
__all__ += functions.__all__
__all__ += precision.__all__
__all__ += ufuncs.__all__

del functions
del precision
del ufuncs