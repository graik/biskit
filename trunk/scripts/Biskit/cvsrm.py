#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

import os
import sys

if (len(sys.argv) < 2):
    print "remove files from disk and cvs"
    print "Syntax: cvsrm file1 file2 file3..."
    print "Changes still need to be committed (cvs ci)."
    sys.exit(0)

for f in sys.argv[1:]:
    os.system("rm "+f)
    os.system("cvs rm "+f)

