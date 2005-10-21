##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$

from os.path import dirname

import tools as T
from Biskit import Trajectory, EHandler
from Biskit.EnsembleTraj import traj2ensemble, EnsembleTraj

from PVM.dispatcher import JobSlave

class QualSlave(JobSlave):
    """
    Plot RMSD, Energy of ensemble trajectory.
    """
    
    def initialize(self, params):

        pass

    def go( self, dict ):

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

        t = T.Load( T.absfile( ftraj ) )
        return traj2ensemble( t )


    def calcRmsd( self, t ):

        mCA = t.ref.maskCA()

        t.fit( ref=t.ref, prof='rms_all_ref', comment='all heavy')

        t.fitMembers( mask=mCA, prof='rms_CA_av', comment='CA to member avg')

        t.fitMembers( refIndex=-1, mask=mCA, prof='rms_CA_last',
                      comment='all CA to last member frame')


    def plotRmsdRef( self, t, title ):
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
