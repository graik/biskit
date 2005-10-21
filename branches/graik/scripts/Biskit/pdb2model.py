#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster
from Biskit.PVM.dispatcher import JobMaster
import Biskit.tools as T
import Biskit.hosts as hosts
import sys

def _use( o ):

    print """
Syntax: pdbs2struct.py -i |file1 file2 ..| [-h |host| -c |chunk| -a -w
                       -o |other_outFolder| -wat -s]

    pvm must be running on the local machine!

Result: pickled PDBModel object for each pdb file with same file name
        but ending in '.model'
        
Options:
    -h    number of hosts to be used
    -a    first add hosts to pvm
    -c    chunk size, number of pdb's passed to each node at once
    -w    display a xterm window for each node
    -o    destination folder (default: same where pdb file comes from)
    -wat  skip water residues (WAT TIP3 WWW H2O)
    -amber  rename CYX -> CYS, HID/HIE/HIP -> HIS, unwrap atom names
          (this creates models with the same atom/res names as pdbs created
           with ambpdb -p top.parm -aatm -bres < some.crd > some.pdb )
    -s    sort atoms alphabetically within residues

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


class StructMaster(TrackingJobMaster):

    def __init__(self, dat, chunk, hosts, outFolder, skipWat=0, amber=0,
                 sort=0, add_hosts=1, **kw):
        """
        dat - data dictionary
        chunk - chunk size
        hosts - list of host names
        outFolder - string, alternative output folder
        """

        niceness = {'default': 0}
        slave_script = T.projectRoot() + '/Biskit/StructureSlave.py'

        TrackingJobMaster.__init__(self, dat, chunk, hosts, niceness,
                           slave_script, **kw)

        self.options = {}
        self.options['out'] = outFolder

        self.options['skipRes'] = None
        if skipWat:
            self.options['skipRes'] = ['WAT','TIP3','H2O','WWW','Na+','Cl-']
        
        if kw.has_key('show_output'):
            self.options['report'] = not kw['show_output']

        self.options['amber'] = amber

        self.options['sort'] = sort


    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.
        """
        return self.options


    def done(self):

        self.exit()
        

if __name__ == '__main__':

    options = T.cmdDict( {'h':10, 'c':5 } )

    if len( sys.argv ) < 2:
        _use( options )

    hostNumber = int( options['h'] )
    chunk = int( options['c'] )
    windows = options.has_key('w')
    amber = options.has_key('amber')
    sort = options.has_key('s')

    outFolder = None
    if options.has_key('o'):
        outFolder = T.absfile( options['o'] )

    ## get all PDBs and models directly from a directory (avoids shell limits)
    if 'd' in options:
        d = T.absfile( options['d'] )
        files = os.listdir( d )
        files = [ x for x in l if x[-4:].upper() == '.PDB' or \
                  x[-6:] == '.MODEL' or x[-7:].upper() == '.PDB.GZ' ]
        files = [ os.path.join( d, f ) for f in files ]
    else:
    ## extract file names
        files = T.toList( options['i'] )

    ## add needed amount of hosts 
    add_hosts = options.has_key('a')

    skipWat = options.has_key('wat')

    fileDic = {}
    for f in files:
        fileDic[T.absfile( f )] = ''

    print "Initialize Job queue.."

    master = StructMaster(fileDic, chunk , hosts.cpus_all[ : hostNumber ] ,
                          outFolder, skipWat=skipWat, amber=amber,
                          sort=sort, show_output=windows, add_hosts=add_hosts )

    print "Start jobs .."
    master.start()

