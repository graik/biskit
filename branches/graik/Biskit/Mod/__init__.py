##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
import Biskit as B

try:

    from Aligner import Aligner
    from Modeller import Modeller
    from SequenceSearcher import SequenceSearcher
    from TemplateSearcher import TemplateSearcher
    from TemplateCleaner import TemplateCleaner

except ImportError, why:
    B.EHandler.warning( 'Error importing Biskit/Mod modules', trace=1 )
