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
"""
Protein-protein docking related modules
"""
from Biskit import EHandler

try:
    from Complex import Complex
    from ComplexEvolving import ComplexEvolving
    from ComplexEvolvingList import ComplexEvolvingList
    from ComplexList import ComplexList
    from ComplexTraj import ComplexTraj
    from ComplexModelRegistry import ComplexModelRegistry
    from ComplexRandomizer import ComplexRandomizer
    from Docker import Docker
    from FixedList import FixedList
    from HexParser import HexParser
    from delphiBindingEnergy import DelphiBindingEnergy
##     from Intervor import Intervor
##     from PatchGenerator import PatchGenerator
##     from PatchGeneratorFromOrbit import PatchGeneratorFromOrbit

except IOError, why:
    EHandler.warning("Couldn't import all Biskit.Dock modules.\n" + str(why))

## PVM-dependent modules
try:
    from ContactMaster import ContactMaster
    from ContactSlave import ContactSlave
except Exception, why:
    EHandler.warning("Couldn't import PVM-dependent modules of Biskit.Dock.\n"+\
                     str( why ) )
