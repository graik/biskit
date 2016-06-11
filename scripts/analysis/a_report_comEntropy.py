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

from Numeric import *
import sys

import biggles as B
import Biskit.tools as T
import Biskit.mathUtils as MU
import Biskit.molUtils as mol

def _use():
    print """
Analyze the result of a_comEntropy.py

Syntax:  a_report_comEntropy.py -i |result.dic| [-eps |file| -tall -tsd -t
                                                 -prefix |str|]
Options:
      -i      dictionary from a_comEntropy.py
      -eps    output plot with entropies for all values of variable parameter
      -tall   print table with all entropy values
      -tsd    print single line with entropy values and standard dev
              (var=='ex3')
      -t      print header row of table
      -prefix prefix for tables (e.g. a01) 

Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def defOptions():
    o = {}
##     o['i'] = '~/interfaces/a11/analysis/entropy/last300_vex3_z500.dat'
    o['o'] = 'test.out'
    o['prefix'] = ''
    return o

class Reporter:

    def __init__( self, tree, protocols ):
        self.tree = tree
        self.prot = protocols
        self.var  = tree['var']
        self.vrange = tree['vrange']

    def arrayDic( self, d, key ):
        r = {}
        for k in d:
            r[k] = array( d[k][key] )
        return r

    def calculate( self, member ):
        """
        d - dict with results of each protocol (for one member or Ensemble)
        -> dict with derrived values
        """
        r = {}
        D = self.arrayDic( self.tree[member], 'S_vibes' )

        ## copy all absolute entropies
        for k in D:
            r['S_' + k] = D[k]
        
        ## binding entropy: com - free 'complex'
        r['dS'] = D['com'] - D['fcom']

        ## total entropy of rec + lig, seperately calculated
        r['S_f'] = D['frec'] + D['flig']
        r['S_b'] = D['brec'] + D['blig']

        ## entropy change internal to rec and lig (i.e. w/o cross-correlations)
        r['dS_intra'] = r['S_b'] - r['S_f']
        r['dS_intra_rec'] = D['brec'] - D['frec']
        r['dS_intra_lig'] = D['blig'] - D['flig']

        ## fake correlations, total
        r['dS_fake1_f'] = D['fcom'] - r['S_f']
        r['dS_fake1_b'] = D['com_split_shift'] - r['S_b']
        r['ddS_fake1']  = r['dS_fake1_b'] - r['dS_fake1_f']

        ## fake correlations, destroyed by shuffling
##         r['dS_fake2_f'] = D['fcom'] - D['fcom_shuff']
##         r['dS_fake2_b'] = D['com_split_shift'] - D['com_split_shuff']
##         r['ddS_fake2'] =  r['dS_fake2_b'] - r['dS_fake2_f']

        ## fake correlations, remaining after shuffling -> MD-stepping
##         r['dS_fake3_f'] = D['fcom_shuff'] - r['S_f']
##         r['dS_fake3_b'] = D['com_split_shuff'] - r['S_b']
##         r['ddS_fake3']  = r['dS_fake3_b'] - r['dS_fake3_f']

        ## binding entropy contributions
        r['dS_rigid'] = D['com'] - D['com_split']             ## rec~lig
        r['dS_cross'] = D['com_split'] - D['com_split_shift'] ## rec#lig
        r['dS_inter'] = D['com'] - D['com_split_shift']    ## rec~lig + rec#lig

        ## entropy calculated without creating fake free complex
        ## add changes of intra and inter-molecular vibrations; correct for
        ## fake correlations arising from regular MD stepping
        r['dS_sep'] = r['dS_inter'] + r['dS_intra'] + r['dS_cross']

        ## translational / rotational entropy loss
        lig_rt = array( self.tree[member]['flig']['S_rot'] ) +\
                 array( self.tree[member]['flig']['S_trans'] )
        rec_rt = array( self.tree[member]['frec']['S_rot'] ) +\
                 array( self.tree[member]['frec']['S_trans'] )
        com_rt = array( self.tree[member]['com']['S_rot'] ) +\
                 array( self.tree[member]['com']['S_trans'] )
        r['S_rt_b'] = com_rt
        r['S_rt_f'] = rec_rt + lig_rt
        r['dS_rt'] = com_rt - ( rec_rt + lig_rt )

        return r

    def table_header( self ):
        s = 15*' ' + 'dS_vibrational\n'
        s += 5*' ' + 32*'-'+'\n'
        s +="%4s %4s  %5s %4s %4s %4s\n" % \
             ('dSrt','total','r|l','r~l','r#l','fake')
        return s

    def table_ensemble( self, member=None, vpos=None, title=1 ):
        """
        Tabulate different entropy components versus variable parameter.
        d      - dict, result of calculate()
        vpos   - int, report result for |i|th variable parameter value
        title  - 1|0, add title lines
        """
        d = self.calculate( member )
        
        s = ''
        if title:
            s += self.table_header()

        vpos = T.toIntList( vpos or range( len( d['dS'] ) ) )

        for i in vpos:
            v = [ d[x][i] for x in ['dS_rt','dS', 'dS_intra', 'dS_rigid',
                                    'dS_cross', 'ddS_fake1'] ]

            s += "%4.0f  %4.0f  %5.0f %4.0f %4.0f %4.0f\n" % \
                 tuple( v )

        return s

    def table_header_error( self, prefix='' ):
        s =  ' '*4 + 23*' ' + 'dS_vibrational' + 23*' ' + 'rottrans\n'
        s += ' '*4 + 58*'-'+'\n'
        s +="    %8s %9s %9s %9s %9s %9s %9s\n" % \
             ('complete','rec', 'lig', 'rec~lig','rec#lig',
              'fake', '')
        return s

    def combined_sd( self, v1, v2 ):
        """
        Calculate the overall standard deviation of two measurements that
        are connected by addition or substraction.
        v1 - [ float ], measurements of value 1
        v2 - [ float ], measurements of value 2
        -> float, standard dev of (v1 +/- v2)
        """
        sd1 = MU.SD( v1 )
        sd2 = MU.SD( v2 )
        return sqrt( sd1**2 + sd2**2 )


    def mean_err( self, d, k_b, k_f ):
        """
        Calculate mean and standard deviation of bound - free entropy from set
        of repeated entropy calculations (leaving out trajectories tripples).
        Mean is in reality the entropy difference using the complete covar.
        matrix. Note that this is a rather arbitrary measure of stability.
        d    - {}, entropy result dictionary
        k_b  - str, key pointing to bound values (i.e. 'com')
        k_f  - str, key pointing to free value (i.e. 'free')
        -> 'mean', 'SD'
        """
        ## S_bound - S_free for all members
        v = d[k_b][0] - d[k_f][0]

        ## get S for remaining variable parameters
        v_f = [ d[ k_f ][i] for i in range(1,len( self.tree['vrange'] ) ) ]
        v_b = [ d[ k_b ][i] for i in range(1,len( self.tree['vrange'] ) ) ]

        return ( v, self.combined_sd( v_f, v_b ) ) 


    def jackknife_diff( self, d, k_b, k_f ):
        """
        Calculate jackknife error estimate of entropy from repeated entropy
        calculations each leaving one trajectory out.
        Note, that jackknife is strictly only defined for the calculation of
        a mean from several samples but is silently applied to all sort of
        other problems.
        d    - {}, entropy result dictionary
        k_b  - str, key pointing to bound values (i.e. 'com')
        k_f  - str, key pointing to free value (i.e. 'free')
        -> full value, error estimate
        """
        ## S_bound - S_free for all members
        v = d[k_b][0] - d[k_f][0]

        ## extract "resampled" entropy differences
        v_resampled = d[ k_b ][1:] - d[ k_f ][1:]

        return v, self.jackknife( v, v_resampled )


    def jackknife( self, v, resampled):
        """
        v         - float, value calculated from complete sample
        resampled - array of float, 'leave one out' resampled values
        -> float, jackknife standard error estimate for v
        """
        n = len( resampled )
        err = sqrt( ((n - 1.0)/n) * sum( (resampled - v)**2 ) )

        return err
        

    def table_ensemble_error( self, title=1, J=0, prefix='' ):
        """
        Use values calculated from variable parameter to determine
        standard dev of first entropy-value.
        """
        d = self.calculate( None )
        
        s = ''
        if title:
            s += self.table_header_error( prefix=prefix )

        x = {}
        x['dS_tot']       = r.jackknife_diff( d, 'S_com', 'S_fcom' )
        x['dS_rt']        = r.jackknife_diff( d, 'S_rt_b', 'S_rt_f')
        x['dS_intra']     = r.jackknife_diff( d, 'S_b', 'S_f' )
        x['dS_rec']       = r.jackknife_diff( d, 'S_brec', 'S_frec' )
        x['dS_lig']       = r.jackknife_diff( d, 'S_blig', 'S_flig' )
        x['dS_rigid']     = d['dS_rigid'][0],\
                            self.jackknife(d['dS_rigid'][0],d['dS_rigid'][1:] )
        x['dS_cross']     = d['dS_cross'][0], \
                            self.jackknife(d['dS_cross'][0],d['dS_cross'][1:] )
        x['dS_fake']      = r.jackknife_diff( d, 'dS_fake1_b', 'dS_fake1_f' )

        ## optionally present dS calculated without fcom..usually worse
        x['dS_sep']       = d['dS_sep'][0], x['dS_intra'][1]

        s += '%3s ' % prefix

        for k in [ 'dS_tot', 'dS_rec', 'dS_lig',
                   'dS_rigid', 'dS_cross', 'dS_fake', 'dS_rt' ]:

            if J==1:
                x[k] = ( mol.calorie * x[k][0], mol.calorie * x[k][1] )
            
            if abs( x[k][0] ) < 10:
                s += "%4.1f" % x[k][0]
            else:
                s += "%4.0f" % x[k][0]

            s += " /%2.0f  " % x[k][1]
                
        s += '\n'

        return s
    

    def plot_ens( self, keys, member=None, firstBold=0, yrange=None ):
        """
        Plot the different entropy components versus variable parameter.
        keys   - [str], keys from R.calculate()
        member - int OR None, member index or None for ensemble
        """
        colors = [ 'black', 'red', 'blue', 'green', 'cyan', 'grey54', 'brown']
        colors.reverse()

        ## axis dimensions and labels
        p = B.FramedPlot()
        p.xlabel = self.var
        p.ylabel = 'bound-free entropy [cal/(mol K)]'

        if yrange: p.yrange = yrange

        d = self.calculate( member )

        sz = firstBold * 3 or 1
        curves = []
        points = []

        for k in keys:
            c = colors.pop()
            curve = B.Curve( self.vrange, d[ k ], color=c, width=sz)
            curve.label = k
            pts = B.Points(self.vrange, d[ k ], color=c, type='filled circle')

            curves.append( curve )
            points.append( pts )

            sz = 1

        ## add first item last (to have it on top)
        curves.reverse(); points.reverse()
        for i in range( len(curves)):
            p.add( curves[i] )
            p.add( points[i] )

        ## Legend
        curves.reverse()
        p.add( B.PlotKey( .70, .85, curves ) )
            
        return p


###############
## MAIN
###############

if __name__ == '__main__':

    options = T.cmdDict( defOptions() )

    if not 'i' in options:
        _use()

    ## Load result tree and dict with protocols
    dat = T.load( options['i'] )
    title = 't' in options

    p = dat['protocols']

    r = Reporter( dat, p )


    if 'eps' in options:
        p = r.plot_ens( ['dS', 'dS_sep', 'dS_rigid', 'dS_cross',
                         'ddS_fake1',
                         'dS_rt'],
                        firstBold=1 )#, yrange=(-200,200) )
        p.show()
        p.write_eps( T.absfile( options['eps'] ) )

    if 'tall' in options:
        print r.table_ensemble()

    if 'tsd' in options:
        print r.table_ensemble_error( J=0, prefix=options['prefix'],
                                      title=title )
