"""
adapted from numpy.oldnumeric version 1.8.0 (dropped from numpy after that). 
"""


from numpy import *

from .functions import *
from .precision import *
from .ufuncs import *
from .misc import *
from .compat import *

from . import functions
from . import precision
from . import ufuncs
from . import misc
from . import compat

__all__ = []
__all__ += functions.__all__
__all__ += precision.__all__
__all__ += ufuncs.__all__
__all__ += misc.__all__
__all__ += compat.__all__

del functions
del precision
del ufuncs
del misc
del compat
