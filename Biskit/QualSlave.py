##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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
## $Revision$
## last $Date$
## last $Author$

"""
Plot RMSD, Energy of ensemble trajectory.
"""
    
from os.path import dirname

import tools as T
from Biskit import Trajectory, EHandler
from Biskit.EnsembleTraj import traj2ensemble, EnsembleTraj

from PVM import JobSlave

class QualSlave(JobSlave):
    """
    Plot RMSD, Energy of ensemble trajectory.
    """

    def initialize(self, params):
        """
        Copy the parameters that Master is passing in as dict into
        fields of this class.
        
        @param params: defined in Master
        @type  params: dict
        """
        pass


    def go( self, dict ):
        """
        Calculate rmsd values for trajectory and plot them.

        @param dict: dictionary with path to trajectories as values
        @type  dict: dict

        @return: dictionary with path to trajectories as values
        @rtype: dict
        """
        for k, f in dict.items():

            fout = dirname( f )+'/'+'%s'+ T.stripFilename(f)+'.eps'
            fshort = dirname( f )[-15:]

            print fshort

            t = self.loadTraj( f )

            self.calcRmsd( t )

            p = self.plotRmsdRef(t, fshort )
            p.write_eps( fout % 'rms_', width="18cm", height="29cm" )

            print "Done"

        return dict


    def loadTraj( self, ftraj ):
        """
        Load trajectories from disc.
        
        @param ftraj: path to trajectory
        @type  ftraj: str

        @return: ensemble trajectory object
        @rtype: EnsembleTraj
        """
        t = T.Load( T.absfile( ftraj ) )
        return traj2ensemble( t )


    def calcRmsd( self, t ):
        """
        Calculate the rmsd to the reference, the CA rmsd to the average
        member and the CA rmsd to last member frame. Add the results
        into a profile.
        
        @param t: ensemble trajectory object
        @type  t:  EnsembleTraj
        """
        mCA = t.ref.maskCA()

        t.fit( ref=t.ref, prof='rms_all_ref', comment='all heavy')

        t.fitMembers( mask=mCA, prof='rms_CA_av', comment='CA to member avg')

        t.fitMembers( refIndex=-1, mask=mCA, prof='rms_CA_last',
                      comment='all CA to last member frame')


    def plotRmsdRef( self, t, title ):
        """
        Plot the rmsd profiles calculated in L{ calcRmsd }.
        
        @param t: ensemble trajectory object
        @type  t:  EnsembleTraj
        @param title: plot title
        @type  title: str

        @return: biggles plot object
        @rtype: biggles.FramedPlot
        """
        p = t.plotMemberProfiles('rms_all_ref', 'rms_CA_av', 'rms_CA_last')
        p.title = title
        p.xlabel= 'frame #'
        p.ylabel= 'RMSD [A]'
        return p


if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = QualSlave()
    slave.start()
