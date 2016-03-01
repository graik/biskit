## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

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
## $Revision$
## last $Date$
## last $Author$
"""
Parallellized AmberEntropist calculation.
"""

import os.path, copy
import numpy as N

import Biskit.tools as T
import Biskit.settings as settings
import Biskit.mathUtils as MU

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster
from Biskit.PVM.hosts import cpus_all, nice_dic
from Biskit import PDBModel, PDBProfiles, EHandler, StdLog

from Biskit.Dock import Complex

slave_path = T.projectRoot()+"/Biskit/AmberEntropySlave.py"

class AmberEntropyMaster(TrackingJobMaster):
    """
    Run many AmberEntropist calculations on many nodes. The Master has
    a standard set of 13 protocols to run on rec, lig, and com
    trajectories, as well as on every single member trajectory - in
    total 113.  It accepts one variable parameter, e.g. s(tart). Each
    protocol is then run for all values of the variable parameter. A
    protocol is simply a set of options that are passed on to the
    AmberEntropist (which is run from within AmberEntropySlave).
    Comparing the different protocols allows to more or less separate
    random from real correlations, rigid body from intermolecular
    vibrations, etc.

    Results are put into a tree-shaped dictionary of dictionaries. The
    first dimension/key is the member index -- None for the complete
    ensemble trajectory, 0 for the first member, etc. The second
    dimension/key is the name of the protocol, e.g. 'com_split' for
    the complex trajectory with seperately fitted receptor and
    ligand. The last dimension contains the different values obtained
    from the ptraj run, e.g. 'S_total' points to the total entropy in
    cal/mol/K, 'contributions' to the entropy contribution of each
    mode, 'T' to the assumed temperature, 'vibes' gives the number of
    vibrations with too low frequencies (according to ptraj). All these
    are lists of values - one for each value of the variable option.

    Example::
             * r[None]['fcom']['S_vibes'][0] -> float
               first vibr. Entropy of free fake complex for complete ensemble
             * r[0]['com']['S_total'] -> [ float, float, .. ]
               the total entropies of the complex calculated for the first
               ensemble member and the different values of the variable option
    """

    def __init__(self, rec=None, lig=None, com=None, out=None,
                 cr=None, var='s', vrange=[0], jack=0,
                 zfilter=None, clean=0, all=1,
                 exrec=[], exlig=[], excom=[],
                 hosts=cpus_all,
                 niceness=nice_dic,
                 w=0, a=1, debug=0,
                 restart=0,
                 **kw ):
        """
        @param rec: free rec trajectory              [required]
        @type  rec: str
        @param lig: free lig trajectory              [required]
        @type  lig: str
        @param com: complex trajectory               [required]
        @type  com: str
        @param out: file name for pickled result     [required]
        @type  out: str
        @param cr: chains of receptor in complex trajectory [n_chains rec]
        @type  cr: [int]

        @param var: name of variable option [ s ]
        @type  var: str
        @param vrange: set of values used for variable option
                       OR 'start:stop:step', string convertable to
                       range() input
        @type  vrange: [any]
        @param jack: set up leave-one-trajectory-out jackknife test
                     (default: 0) (replaces var with 'ex1' and vrange with
                     range(1,n_members+1))
        @type  jack: [0|1]

        @param zfilter: kick out outlyer trajectories using z-score threshold
                        on RMSD trace (default: None->don't)
        @type  zfilter: float
        @param clean: remove pickled ref models and member trajectories
                      (default: 0)
        @type  clean: 0|1
        @param all: skip single member trajs (default: 1)
        @type  all: 0|1

        @param exrec: exclude certain members of receptor ensemble    [[]]
        @type  exrec: [int]
        @param exlig: exclude certain members of ligand  ensemble     [[]]
        @type  exlig: [int]
        @param excom: exclude certain members of complex ensemble     [[]]
        @type  excom: [int]

        @param hosts: nodes to be used (default: all known)
        @type  hosts: [str]
        @param debug: don't delete output files (default: 0)
        @type  debug: 1|0


        @param kw: additional key=value parameters for AmberEntropist,
                   AmberCrdEntropist, Executor and Master.
        @type  kw: key=value pairs
        ::
          ... parameters for AmberEntropist
          cast    - 1|0, equalize free and bound atom content [1]
          s,e     - int, start and stop frame                 [0, to end]
          atoms   - [ str ], names of atoms to consider       [all]
          protein - 1|0, remove non-protein atoms             [0..don't]
          step    - int, frame offset                         [no offset]
          thin    - float, use randomly distributed fraction of frames [all]
                    (similar to step but perhaps better for entropy
                    calculations)
          ex      - [int] OR ([int],[int]), exclude member trajectories   [[]]
          ex_n    - int, exclude last n members  OR...                  [None]
          ex3     - int, exclude |ex3|rd tripple of trajectories          [0]
                    (index starts with 1! 0 to exclude nothing)

          ... parameters for AmberCrdEntropist
          f_template - str, alternative ptraj input template  [default]
  
          ... parameters for Executor:
          log      - Biskit.LogFile, program log (None->STOUT)        [None]
          verbose  - 0|1, print progress messages to log     [log != STDOUT]
 
          ... parameters for Master
          w        - 0|1, show X window for each slave [0]
          a        - 0|1, add hosts to PVM [1]
        """
        ## normal and error output
        self.fout = T.absfile( out )
        self.ferror = os.path.dirname(self.fout) +'/AmberEntropy_errors.log'
        self.debug = debug

        self.log = StdLog()

        ## input files and variable option
        self.rec = T.absfile( rec, 0 )
        self.lig = T.absfile( lig, 0 )
        self.com = T.absfile( com, 0 )
        self.cr  = cr
        self.cl  = None
        self.var    = var
        self.vrange = self.__vrange( vrange )
        self.jack = jack
        self.zfilter = zfilter
        self.n_members = None
        self.clean = clean
        self.all = all

        ## members to exclude, outliers will be added if zfilter != None
        self.ex_frec = exrec
        self.ex_flig = exlig
        self.ex_com  = excom

        ## reserve for loaded reference models
        self.ref_frec = self.ref_flig = None
        self.ref_brec = self.ref_blig = self.ref_com  = None

        ## reserve for extracted member trajectories
        self.members_frec = self.members_flig = []
        self.members_brec = self.members_blig = []

        ## options to be passed on to AmberEntropist
        self.options = kw

        if not restart:
            ## Load trajectories, dump references, identify outliers
            self.processTrajs()

            ## prepare dict of protocols for AmberEntropist
            self.protocols = self.protocols_var_range( **kw )
            self.saveProtocols()

        TrackingJobMaster.__init__(self, self.protocols,
                                   chunk_size=1,
                                   hosts=hosts,
                                   niceness=niceness,
                                   slave_script=slave_path,
                                   show_output=w,
                                   add_hosts=a)

        print "JobMaster initialized."


    def __vrange( self, v ):
        """
        Interprete the vrange option -> [ int ] or [ float ]

        @param v: vrange option
        @type  v: lst OR str
        
        @return: range option
        @rtype: [int] OR [float]
        """
        if type( v ) is list:
            return [ self.__float_int(x) for x in v ]
        if type( v ) is str and ':' in v:
            v = tuple( [ self.__float_int(x) for x in v.split(':') ] )
            return N.arange( *v )

        return self.__float_int( v )


    def __float_int( self, v ):
        """
        Convert v to int or, if necessary, float

        @param v: value
        @type  v: any

        @return: converted value
        @rtype: int OR float        
        """
        if float(v) % 1. != 0:
            return float( v )
        return int( float(v) )


    def loadTraj( self, fname, outliers=[], refname=None  ):
        """
        Load trajectory from file.
        
        @param fname: path to trajectory
        @type  fname: str
        @param outliers: Identify outlier trajectories (default: [], identify)
        @type  outliers: [int] OR [] 
        @param refname: name of reference (efault: None)
        @type  refname: str
        
        @return: t, outliers, members
        @rtype: trajectoty, [int], [int]
        """
        self.log.add('Loading ' + fname )
        t = T.load( fname )

        t.ref.addChainId()
        t = t.compressAtoms( t.ref.maskProtein() )

        outliers = self.getOutliers( t, outliers )

        if refname:
            self.dumpMissing( t.ref, refname )

        members = None
        if not self.all:
            members = self.dumpMembers( t, self.rec )

        return t, outliers, members


    def processTrajs( self ):
        """
        Extract reference model and member trajectories from rec, lig, and
        com trajectories. Identify outlier member trajectories, if requested.
        """
        ## free rec
        self.ref_frec = self.nameRef( self.rec )

        t, self.ex_frec, self.members_frec = self.loadTraj(
            self.rec, self.ex_frec, self.ref_frec )

        n_rec_members = t.n_members
        self.cr       = self.cr or range( t.ref.lenChains( breaks=0 ) )
        del t

        ## free lig
        self.ref_flig = self.nameRef( self.lig )

        t, self.ex_flig, self.members_flig = self.loadTraj(
            self.lig, self.ex_flig, self.ref_flig )

        n_lig_members = t.n_members
        del t

        ## complex
        fname = T.stripSuffix( T.absfile( self.com, resolveLinks=0 ) )
        self.ref_com = fname + '_ref.complex'
        self.ref_blig= fname + '_blig.model'
        self.ref_brec= fname + '_brec.model'

        t, self.ex_com, self.members_com = self.loadTraj(
            self.com, self.ex_com  )

        n_com_members = t.n_members

        self.cl = self.cl or MU.difference( range(t.ref.lenChains()), self.cr)
        rec = t.ref.takeChains( self.cr, breaks=0 )
        lig = t.ref.takeChains( self.cl, breaks=0 )

        del t
        self.dumpMissing( Complex( rec, lig ), self.ref_com )
        self.dumpMissing( rec, self.ref_brec )
        self.dumpMissing( lig, self.ref_blig )

        self.equalizeMemberCount( n_rec_members, n_lig_members, n_com_members )

        if self.jack:  self.prepareJackknife()


    def equalizeMemberCount( self, n_rec, n_lig, n_com ):
        """
        Ensure we keep equal number of members trajectories from frec,
        flig, and com.

        @param n_rec: number of receptor trajectories
        @type  n_rec: int
        @param n_lig: number of ligand trajectories
        @type  n_lig: int
        @param n_com: number of complex trajectories
        @type  n_com: int
        """
        ex        = [ self.ex_frec, self.ex_flig, self.ex_com ]
        n_members = [ n_rec, n_lig, n_com ]

        ## pair list of excluded members with number of remaining members
        ex = [ ( ex[i], n_members[i] - len(ex[i]) ) for i in range(3) ]

        ## lowest number of members after exclusion
        n_min = min( [ x[1] for x in ex ] )

        self.log.add('excluding non-outliers to match member count: ')

        label = ['com','lig','rec']

        for x, n in ex:
            i = 0
            s = label.pop()

            while n > n_min:
                self.log.write( '%s: ' % s )
                if not i in x:
                    x.append( i )
                    n -= 1
                    self.log.write('%i, ' % i )
                i += 1

            self.log.add('')

        self.n_members = n_min


    def prepareJackknife( self ):
        """
        Prepare leave-one-trajectory-out jackknife test.
        """
        self.vrange = range( self.n_members + 1 )  ## 0: exclude nothing
        self.var = 'ex1'


    def nameRef( self, fname ):
        fname = T.stripSuffix( T.absfile( fname, resolveLinks=0 ) )
        return fname + '_ref.model'


    def nameRefCom( self, fname ):
        fname = T.stripSuffix( T.absfile( fname, resolveLinks=0 ) )
        return fname + '_ref.complex'


    def dumpMissing( self, o, fname ):
        """
        Pickle *o* to path *fname*, if it is not already there.

        @param o: object to dump
        @type  o: any
        @param fname: file name
        @type  fname: str
        
        @return: file name
        @rtype: str
        """
        if os.path.exists( fname ):
            self.log.add('using existing ' + fname )
        else:
            self.log.add('Saving ' + fname )
            T.dump( o, fname )

        return fname


    def getOutliers( self, traj, outlaws=[] ):
        """
        Identify member trajectories that haved moved much further than normal.

        @param traj: Trajectory to analyze
        @type  traj: Trajectory
        @param outlaws: members already marked for exclusion
        @type  outlaws: [int]

        @return: member indices of outlyer trajectories (plus outlaws)
        @rtype: [int]
        """
        if not self.zfilter:
            return outlaws

        outliers = N.nonzero( traj.outliers( z=self.zfilter,
                                             mask=traj.ref.maskCA(), step=10) )[0]
        self.log.add('identified %i outliers with z-threshold %3.1f' %\
                     ( len(outliers), self.zfilter ) )

        return MU.union( outliers, outlaws )


    def dumpMembers( self, traj, fname  ):
        """
        Dump ensemble member trajectories
        
        @param traj: Trajectory to dump
        @type  traj: Trajectory
        @param fname: trajectory file name - used to derrive name for members
        @type  fname: str'
        
        @return: list of trajectory files
        @rtype: [str]
        """
        fname = T.stripSuffix( T.absfile( fname, resolveLinks=0 ) )
        members = range( traj.n_members )

        r = []
        for n in members:
            f = fname + '_member_%02i.traj' % n
            if os.path.exists( f ):
                self.log.add('using existing ' + f )
            else:
                self.log.write('saving ' + f + '...')
                m = traj.takeMember( n )
                T.dump( m, f )
                self.log.add('done')
            r += [ f ]

        return r


    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}
        """
        host = self.hostnameFromTID( slave_tid )
        nice = self.niceness.get( host, self.niceness.get('default',0) )

        return {'ferror':self.ferror,
                'debug':self.debug, 'nice':nice, 'host':host}


    def cleanup( self ):
        """
        Tidy up
        """
        if self.clean:
            self.cleanCache()

    def cleanCache( self ):
        """
        Remove left-over cache files
        """
        fs = [ self.ref_frec, self.ref_flig, self.ref_com, self.ref_brec,
               self.ref_blig ]
        fs.extend( self.members_frec + self.members_flig )
        fs.extend( self.members_brec + self.members_blig )
        fs.extend( self.members_com )

        for f in fs:
            self.log.add('removing %s: %i' % (f, T.tryRemove(f)) )


    def saveProtocols( self ):
        """
        Save protocol to file.
        """
        f_prot = T.stripSuffix( T.absfile(self.fout) ) + '_protocols.dat'
        self.log.write( 'Saving parameters to %s...' % f_prot )
        T.dump( self.protocols, f_prot )


    def done(self):
        """
        Write result to file.
        """
        tree = self.getResult()
        self.log.add("Saving result to %s..." % self.fout)
        T.dump( tree, self.fout )
        self.log.add( "Done" )


    ##
    ## Assemble the protocols for many AmberEntropist runs
    ##
    def __cpupdate( self, d1, d2 ):
        """
        Merge 2 dictionaries *d1* and *d2* and return a copy
        """
        r = copy.copy( d1 )
        r.update( d2 )
        return r

    def protocols_standard( self, trec, tlig, tcom,
                            ex_frec=None, ex_flig=None, ex_com=None,
                            doshift=1,
                            **options ):
        """
        Create 13 parameter sets for AmberEntropist that cover the calculation
        of rec, lig, com and fcom entropies with and without splitting of the
        complex, with and without shifting and shuffling of frames.
        
        @param options: additional options (like cast, s, e, atoms, thin, step)
                        that are the same in all parameter sets
        @type  options: key=value
        
        @return: each value of the returned dict contains a set of
                 arguments for one AmberEntropist run
        @rtype: dict of dict
        """
        fcp = self.__cpupdate
        r = {}
        S = self  ## make rest more readable

        d = { 'ref':None, 'cast':1, 'chains':None,
              'split':0, 'shift':0, 'shuffle':0, 'ex_n':0, 'ex3':None,
              'thin':None, 'step':1, 'ss':0, 'se':None, 'atoms':None }
        d.update( options )

        r['frec'] = fcp( d, {'traj':trec, 'ref':S.ref_brec, 'ex':ex_frec } )
        r['flig'] = fcp( d, {'traj':tlig, 'ref':S.ref_blig, 'ex':ex_flig } )
        r['brec'] = fcp( d, {'traj':tcom, 'ref':S.ref_frec, 'ex':ex_com,
                             'chains':S.cr } )
        r['blig'] = fcp( d, {'traj':tcom, 'ref':S.ref_flig, 'ex':ex_com,
                             'chains':S.cl } )

        r['fcom'] = fcp( d, {'traj':'%s+%s'%(trec, tlig),
                             'ex':(ex_frec, ex_flig),
                             'ref':S.ref_com, 'split':1 } )

##         if doshift:
##             r['fcom_shift'] = fcp( r['fcom'], {'shift':1 } )

        r['fcom_shuff'] = fcp( r['fcom'], {'shuffle':1 } )

        r['com']  = fcp( d, {'traj':tcom, 'ex':ex_com,
                             'ref':'%s+%s' % (S.ref_frec, S.ref_flig) } )

        r['com_split'] = fcp( r['com'], { 'split':1,   'border':S.cl[0] } )
##      r['com_shuff'] = fcp( r['com'], { 'shuffle':1, 'border':S.cl[0] } )
        r['com_split_shuff'] = fcp( r['com'],
                                {'split':1,'shuffle':1,'border':S.cl[0] } )
        if doshift:
##             r['com_shift'] = fcp( r['com'], { 'shift':1,'border':S.cl[0] } )
            r['com_split_shift'] = fcp( r['com'],
                                     {'split':1,'shift':1, 'border':S.cl[0] } )

        return r


    def protocols_single_all( self, **options ):
        """
        Set of protocols for all-member trajectories AND single-member traj.
        with the different shuffle, shift, split settings.
        Usually 11 x 13 protocols for AmberEntropist (10 members and 1 for all)
        
        @param options: additional options (like cast, s, e, atoms, thin, step)
                        that are the same in all parameter sets
        @type  options: key=value
        
        @return: each value of the returned dict contains a set of arguments
                 for one AmberEntropist run, each key is a tuple of the
                 member index and the protocol name, i.e. (0, 'fcom_shuffle')
                 The set of protocols for all-member trajectories has member
                 index None.
        @rtype: dict of dict
        """
        r = {}
        ## put all-member protocolls under member index 'None'
        prots = self.protocols_standard( self.rec, self.lig, self.com,
                                     self.ex_frec, self.ex_flig, self.ex_com,
                                     **options )
        for k,p in prots.items():
            r[ (None, k) ] = p 

        if not self.all:
        ## put single-member protocols under their respective member index
            for i in range( len( self.members_frec ) ):
                prots = self.protocols_standard(self.members_frec[i],
                                            self.members_flig[i],
                                            self.members_com[i], doshift=0,
                                            **options )
                for k, p in prots.items():
                    r[ (i, k) ] = p

        return r


    def protocols_var_range( self, **options ):
        """
        Complete set of protocols also considering different values of the
        variable option.
        """
        self.log.add( 'variable option %s with %i values' \
                      % (self.var, len(self.vrange)))

        r = {}
        for v in self.vrange:
            d = copy.copy( options )
            d[ self.var ] = v

            prots = self.protocols_single_all( **d )

            for k, p in prots.items():
                r[ (v,) + k ] = p

        return r

    ##
    ## Re-organize results
    ##
    def dictionate( self, d ):
        """
        Take dict with tuple keys (value, int_member, str_protocol) and build
        a tree-like dict of dicts in which the values of d can be accessed
        like::
          d[value][int_member][str_protocol]

        @param d: the raw results accumulated from the slave nodes
        @type d: dict

        @return: tree-like dict ordered by variable value, member, protocol
        @rtype: dict of dict of dict of dict
        """
        r = {}

        keys = d.keys()

        ## only convert single value tuple keys into non-tuple keys
        if len( keys[0] ) == 1:
            for k in keys:
                r[ k[0] ] = d[ k ]
            return r

        x_values = MU.nonredundant( [ k[0] for k in keys ] )

        for x in x_values:

            sub_keys = [ k for k in keys if k[0] == x ]
            y_values = MU.nonredundant( [ k[1:] for k in sub_keys] )

            r[ x ] = {}
            for y in y_values:
                r[x][y] = d[ (x,) + y ]

            r[ x ] = self.dictionate( r[x] )

        return r


    def getResult( self, **arg ):
        """
        Collapse the results for different values of the variable parameter
        into lists and put the results into a tree ala::
          r[ member_index ][ protocol_name ][ result_field ] -> [ values ]

        @return: tree-like dict ordered by variable value, member, protocol
        @rtype: dict of dict of dict of lists
        """
        tree = self.dictionate( self.result )

        vvalues = tree.keys()
        vvalues.sort()

        keys = self.result.keys()
        sub_keys = [ k for k in keys if k[0] == vvalues[0] ]

        r = {}
        for v, member, protcl in sub_keys:

            try:
                if not member in r:
                    r[member] = {}

                r[member][protcl] = {}

                run_dic = tree[v][member][protcl]

                for k in run_dic.keys():
                    r[member][protcl][k] = [ tree[v][member][protcl][k] \
                                             for v in vvalues ]
            except:
                EHandler.warning('missing result: ' + str(T.lastError()))

        r['var'] = self.var
        r['vrange']= self.vrange
        r['protocols'] = self.protocols

        self.result_tree = r
        return r




#### TEST #######

if __name__ == '__main__':
    niceness = {'default': 0}
    hosts = cpus_all[:80]

    f = T.testRoot() + '/Amber/AmberEntropyMaster/'

    rec = f + 'rec/traj.dat'
    lig = f + 'lig/traj.dat'
    com = f + 'com/traj.dat'
    out = f + 'entropy.out'

    master = AmberEntropyMaster( rec, lig, com, out, step=1,
                                 atoms=['CA','CB'],
                                 var='ex1', vrange='0:10',
                                 exrec=[1],exlig=[0],
                                 all=1,
                                 hosts=hosts, niceness=niceness,
                                 w=1 )

    master.start()
