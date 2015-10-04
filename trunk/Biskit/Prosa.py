## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
## Author:        Wolfgang Rieping
## Created:       04/02/02
## Last modified: 04/06/02
##
"""
Run a ProsaII job.

@attention: This class should be replaced by L{Prosa2003}. It will only stay
for a while unitl  all transitions are made to the new version of Prosa.
"""

import tempfile
import os
import numpy as N

import tools as T
import Table
import settings as S

class ProsaII:

    def __init__(self, executable = 'prosaII', temp_dir=S.tempDirShared ):

        self.executable = executable
        self.setCleanup()
        self.temp_dir = temp_dir
        self.script_name = tempfile.mktemp()

        ## set default values

        self.setPairWindowSize()
        self.setSurfaceRange()
        self.setSurfaceFactor()

    def setCleanup(self, value = 1):

        self.cleanup = int(value)

    def getCleanup(self):

        return self.cleanup

    def setPairWindowSize(self, window = (1, 600)):
        """
        pair-interactions are calculated for residue pairs
        whose sequence separation lies in 'window'.
        window must be a tuple.
        default values come from prosaII manual 
        """

        self.pairWindowSize = tuple(window)

    def getPairWindowSize(self):

        return self.pairWindowSize

    def setSurfaceRange(self, interval = (0., 15.)):
        """
        surface energy for a residue pair is calculated only,
        if their distance lies within the 'interval'. otherwise
        energy is 0.
        interval must be tuple.
        default values come from prosaII manual
        """

        self.surfaceRange = tuple(interval)

    def getSurfaceRange(self):

        return self.surfaceRange

    def setSurfaceFactor(self, factor = 5.):
        """
        E_total = factor * E_surface + E_pair
        see manual for further recommendations
        """

        self.surface_factor = factor

    def getSurfaceFactor(self):

        return self.surface_factor

    def run(self, command):

        script = self.script_name

        f = open(script, 'w')
        f.write(command)
        f.close()

        os.system('%s -d -f %s > /dev/null' % (self.executable, script))

        if self.getCleanup():

            try:
                os.unlink(script)
            except:
                pass

    def analyseEnergy(self, filename, object_name = None):
        """
        supports user-expansion
        """

        if not os.path.exists(filename):
            raise IOError, 'file %s does not exist. Check the ProsaII license' % filename

        filename = os.path.expanduser(filename)
        path, filename = os.path.split(filename)
        name, ext = os.path.splitext(filename)

        if path == '':
            path = '.'

        lower_k, upper_k = self.getPairWindowSize()
        surface_lower, surface_upper = self.getSurfaceRange()

        if object_name is None:
            object_name = name

        values = {'pdb_path': path,
                  'pdb_file': filename,
                  'obj_name': object_name,
                  'lower_k': lower_k,
                  'upper_k': upper_k,
                  'surface_lower': surface_lower,
                  'surface_upper': surface_upper,
                  'factor_surface': self.getSurfaceFactor()}

        command = 'pdb_dir = %(pdb_path)s\n' + \
                  'read pdb %(pdb_file)s %(obj_name)s\n' + \
                  'lower_k = %(lower_k)d\n' + \
                  'upper_k = %(upper_k)d\n' + \
                  'pot_lb = %(surface_lower)f\n' + \
                  'pot_ub = %(surface_upper)f\n' + \
                  'analyse energy %(obj_name)s\n' + \
                  'print energy %(obj_name)s\n'

        ## save current path

        old_path = os.getcwd()
        os.chdir(self.temp_dir)

        self.run(command % values)

        ## output is in .ana-file

        prosa_output = values['obj_name'] + '.ana'

        result = parse_ANA_output(prosa_output)

        ## if cleanup-flag is set, remove ana-files

        if self.getCleanup():
            os.unlink(prosa_output)

        ## restore old path

        os.chdir(old_path)

        return result[:,1:]


def parse_ANA_output(prosa_output):

    if not os.path.exists(prosa_output):
        raise 'prosa output %s is missing.' % prosa_output

    t = Table.fromFile(prosa_output)

    return N.array(t[2:]).astype(N.float32)

#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """no test assigned"""

    TAGS = [ BT.OLD, BT.EXE ]

    pass
