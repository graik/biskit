##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##
import Biskit as B

try:

    from Aligner import Aligner
    from Modeller import Modeller
    from SequenceSearcher import SequenceSearcher
    from TemplateSearcher import TemplateSearcher
    from TemplateCleaner import TemplateCleaner
    from TemplateFilter import TemplateFilter

    import Biskit.PVM as PVM
    if PVM.pvm_installed:
        from ModelMaster import ModelMaster
        from ModelSlave  import ModelSlave
        from AlignerMaster import AlignerMaster
        from AlignerSlave  import AlignerSlave

except ImportError, why:
    B.EHandler.warning( 'Error importing Biskit/Mod modules', trace=1 )

try:
    del PVM
except:
    pass
