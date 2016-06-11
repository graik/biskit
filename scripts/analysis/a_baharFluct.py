from a_baharEntropy import *
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
import biggles

def go( f, fa ):

    ## get lattice model CA fluctuation
    o = {}
    o['com'] = f + '/com_wet/dry_com.model'
    o['rec'] = f + '/rec_wet/dry_rec.model'
    o['lig'] = f + '/rec_wet/dry_rec.model'

    a = BaharAnalyzer( **o )

    a.getMatrices()

    fbahar = a.fluctuation( a.m_b )

    ## get xplor and amber trajectories
    t_cx = load(f  + '/com_pcr_00/traj_fluct.dat')
    t_ca = load(fa + '/com_pcr_00/traj_fluct.dat')

    t_cx = t_cx.compressAtoms( t_cx.ref.maskCA() )
    t_ca = t_ca.compressAtoms( t_ca.ref.maskCA() )

    ## get CA fluct profiles for Xplor and Amber Traj
    fx = t_cx.ref.profile('fluct_global')
    fa = t_ca.ref.profile('fluct_global')

    ## Plotting

    p = biggles.FramedPlot()
    p.add( biggles.Curve( range( len( fx ) ), fx, color='blue' ) )
    p.yrange = (0,1.5)

    p.add( biggles.Curve( range( len( fa ) ), fa, color='black') )

    p.add( biggles.Curve( range( len( fbahar)), fbahar, color='red') )

    return p

########## MAIN ##############

p11 = go( '~/interfaces/c11', '~/interfaces/a11' )
p02 = go( '~/interfaces/c02', '~/interfaces/a02' )
p17 = go( '~/interfaces/c17', '~/interfaces/a17' )

p = biggles.Table( 3, 1 )
p[0,0] = p02
p[1,0] = p11
p[2,0] = p17
