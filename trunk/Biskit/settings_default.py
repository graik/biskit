##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## Default parameters for Biskit.settings
##
## last $Author$
## last $Date$
## $Revision$

"""
Defines the options exported by settings.py and 
provides fall-back values for options that are not found in
.biskit/settings.dat.

This is NOT the place where you adapt biskit to the local environment!
Local settings belong into .biskit/settings.dat which is then parsed by
Biskit/settings.py. 

Nevertheless, this IS the place where you define new global parameters for
biskit. These parameters become then (semi-magically) available as static
fields of Biskit.settings. settings.py is importing the content of this module.
For each field found here, it expects an option in .biskit/settings.dat
under the given section. If such an option (or .biskit/settings.dat) is not
found, the default parameter is taken (or None) and a Warning is issued.

Each parameter is defined by a tuple containing
  (1) the default value for a parameter of this name in settings.py
  (2) the section into which this parameter is put within .biskit/settings.dat
  (3) an optional comment (which also ends up in .biskit/settings.dat)

Note: Keep execution time of this module to a minimum! It always needs to
be imported at the beginning of a Biskit session. That means, stay away
from file searching or similar things.
"""

import os as __os
import user as __user
import Biskit.tools as __T

## cache for speed
__root = __T.projectRoot()

#################################################
## Parameters expected in .biskit/settings.dat ##

## BISKIT_BIN
s = 'biskit_bin' ## current section (for .biskit/settings.dat)

ssh_bin          = '/usr/bin/ssh', s, 'Path to ssh executable'
nice_bin         = '/bin/nice', s
whatif_bin       = 'whatif', s, 'WhatIF: http://swift.cmbi.kun.nl/whatif'
prosa2003_bin    = 'prosa2003', s,\
                   'http://www.proceryon.com/solutions/prosa.html'
fold_x_bin 	 = __root + '/external/fold_x_linux', s,\
                    'FoldX: http://foldx.embl.de'

bl2seq_bin       = 'bl2seq', s,\
                   'part of NCBI-Blast; http://www.ncbi.nlm.nih.gov/BLAST/'

hmmfetch_bin     = 'hmmfetch', s, 'Part of HMMER; http://hmmer.wustl.edu'
hmmpfam_bin      = 'hmmpfam', s,  'Part of HMMER; http://hmmer.wustl.edu'
hmmalign_bin     = 'hmmalign', s, 'Part of HMMER; http://hmmer.wustl.edu'
hmmindex_bin     = 'hmmindex', s, 'Part of HMMER; http://hmmer.wustl.edu'

msms_bin 	 = 'msms', s, 'http://www.scripps.edu/~sanner/python/'
pdb_to_xyzrn_bin = __root + '/external/msms/gpdb_to_xyzrn', s,\
                    'part of the the MSMS package'

fastsurf_bin      = 'fastsurf7_2', s,\
                    'http://monte.biochem.wisc.edu/~tsodikov/surface.html'
xplor_bin 	  = 'ifcxplor', s, 'X-Plor with PCR-restraint capabilities'
gbxplor_bin       = 'gbxplor', s, 'X-Plor with Generalized Born implementation'

pymol_bin 	  = 'pymol.com', s, 'http://pymol.sourceforge.net/'

##icmbrowser_bin    = 'icmbrowser', s, 'http://www.molsoft.com/icm_browser.html'

## AMBER
s = 'amber'

ptraj_bin         = 'ptraj', s,  'part of Amber package'
tleap_bin         = 'tleap', s,  'part of Amber package'
sander_bin        = 'sander', s, 'part of Amber package'
ambpdb_bin        = 'ambpdb', s, 'part of Amber package'
leaprc            = __os.getenv('AMBERHOME','')+'/dat/leap/cmd/leaprc.ff99',s,\
                    'default force field for Amber tleap'


## BISKIT_PATHS
s = 'biskit_paths'
tempDirLocal      = __T.tempDir(), s, 'Local Temp directory (default: /tmp)'
tempDirShared     = __user.home + '/tmp', s,\
                    'Tmp directory reachable by all hosts'
pymol_scripts     = __root + '/external/pymol/', s,\
                    'external scripts for PyMol'


## BISKIT_DB
s = 'biskit_db'

hmm_db            = '/Bis/db/pfam/hmmer/Pfam', s,\
                    'Pfam (http://pfam.wustl.edu) for HMMER'

s = 'test'

test_bin    = 'superprogram',   s, 'test parameter, not used'
test_float  = 1.2, s, 'test parameter, not used'

#############################
## clean up module content ##

del s, __root, __os, __T, __user
