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
        self.verbose=params.get('verbose',1)


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

            if self.verbose:
                print fshort

            t = self.loadTraj( f )

            self.calcRmsd( t )

            p = self.plotRmsdRef(t, fshort )
            p.write_eps( fout % 'rms_', width="18cm", height="29cm" )

            if self.verbose:
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
        t = T.load( T.absfile( ftraj ) )
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

        t.fit( ref=t.ref, prof='rms_all_ref', comment='all heavy',
               verbose=self.verbose )

        t.fitMembers( mask=mCA, prof='rms_CA_av', comment='CA to member avg',
                      verbose=self.verbose )

        t.fitMembers( refIndex=-1, mask=mCA, prof='rms_CA_last',
                      comment='all CA to last member frame',
                      verbose=self.verbose)


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

#############
## DEBUGGING
#############
from threading import Thread
import time

class Flusher( Thread ):

    def __init__( self, *f ):
        Thread.__init__( self )
        self.setDaemon( True )

        self.files = f
        self.stop = False

    def setStop( self ):
        self.stop = True

    def run( self ):

        while not self.stop:
            for f in self.files:
                f.flush()
            time.sleep( 2 )


def redirect_output():
    import tempfile, os, sys

    f_out = open( tempfile.mktemp( '.out', 'slave_', T.absfile('~')),'w')

    sys.stdout = f_out
    sys.stderr = f_out
    flusher = Flusher( f_out )
    flusher.start()



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test QualSlave without the master."""

    TAGS = [ BT.PVM ]

    def test_QualSlave(self):
        """QualSlave test"""
        import os

        jobs = {0: T.testRoot() + '/lig_pcr_00/traj.dat'}

        self.feps = '%s/rms_traj.eps' % os.path.dirname( jobs[0] )

        slave = QualSlave()
        slave.initialize( {'verbose':0} )
        r = slave.go( jobs )

        if self.VERBOSITY > 2:
            self.log.write('eps written to ' + self.feps )

        self.assert_( os.path.exists( self.feps ) )

    def cleanUp(self):
        T.tryRemove( self.feps )

if __name__ == '__main__':

##     BT.localTest()  ## for interactive debugging only

    redirect_output()

    import os, sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = QualSlave()
    slave.start()
