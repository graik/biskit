##
## Biskit, a toolkit for the manipulation of macromolecular structures
## coiledcoil -- a package for coiled-coil analysis
## Copyright (C) 2009 Victor Gil
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
## $Revision: 652 $
## last $Author: graik $
## last $Date: 2009-03-03 07:38:23 +0100 (Tue, 03 Mar 2009) $

"""
Analysis and alignment of coiled-coil sequences.
"""
from Biskit import EHandler

## import the public API classes into package namespace
try:
    from coiledcoil import CoiledCoil
    from coiledalign import CoiledAlign
except IOError, why:
    EHandler.warning("Couldn't import all Biskit.coiledcoil modules.\n" + str(why))
