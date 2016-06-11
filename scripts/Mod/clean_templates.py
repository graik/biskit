#!/usr/bin/env python
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
##

import Biskit.Mod.modUtils as modUtils
from Biskit.Mod import *
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile

import sys, os.path

def _use( o ):
    print """
Syntax: clean_templates.py [-o |output_folder| -i |chainIndex| -log |logFile|

input: templates/nr/*.pdb
       templates/nr/chain_index.txt

output: templates/t_coffee/*.alpha    (input for Alignar)
        templates/modeller/*.pdb      (input for Modeller)
        
Options:
    -o       output folder for results      (default: .)
    -i       chain index file for templates
                 (default: '/templates/nr/chain_index.txt')
    -log     log file                       (default: STDOUT)

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value
    sys.exit(0)


def defaultOptions():
    return {'o':'.',
            'log': None
            }

### MAIN ###

options   = tools.cmdDict( defaultOptions() )
outFolder = tools.absfile( options['o'] )
chIndex = options.get('i', outFolder + TemplateSearcher.F_NR +
                      TemplateSearcher.F_CHAIN_INDEX )

chIndex = tools.absfile( chIndex )

if not (os.path.exists( chIndex ) ):
    _use( defaultOptions() )

if '?' in options or 'help' in options:
    _use( defaultOptions() )

log = None
if options['log']:
    log = LogFile( outFolder + '/' + options['log'], 'a' ) 

###################
## TemplateCleaner
##
## Prepare pdb files in templates/nr for T-coffee and modeller
## (replace nonstandard residues, remove hydrogens,
#    remove atoms with nultiple configurations, etc.)

## input: templates/nr/*.pdb
##        templates/nr/chain_index.txt
##
## output: templates/t_coffee/*.alpha    (input for Alignar)
##         templates/modeller/*.pdb      (input for Modeller)

try:
    cleaner = TemplateCleaner( outFolder, log )

    inp_dic = modUtils.parse_tabbed_file( chIndex )

    cleaner.process_all( inp_dic )

except:
    EHandler.error( 'Error while cleaning templates')

