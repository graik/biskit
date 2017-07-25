#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg
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

import sys, copy

import Biskit.tools as T
from Biskit.Dock import DelphiBindingEnergy, Complex
from Biskit import PDBModel

def _use( options ):
    print """
delphiBinding.py - calculate the electrostatic component of binding free Energy
                   with the Delphi PB solver. The input complex object can be
                   created with pdb2complex.py. This script is a wrapper for
                   the Biskit.Dock.DelphiBindingEnergy class and requires
                   two external programs: (1) reduce and (2) delphi. 

Syntax:  delphiBinding.py  -c |input complex|
                            OR
                           -i |input PDB file| 
                           -r |chain index| -l |chain index|

                           -o |result file|
                           -ocom |output complex|
                           -log |log file|
                           -v 
                           -debug

                           -autocap
                           -indi |float| -exdi |float| -salt |float|
                           -ionrad |float| -prbrad |float| -bndcon |float| 
                           -scale |float| -perfil |int|
                           
                     
Options:   -c     complex pdb file or pickled Biskit.Dock.Complex object
           -r     receptor chain list (e.g. 0 1 ) ignored if Complex is given
           -l     ligand      ~       (e.g. 2 )   ignored if Complex is given
           -o     output text file with energy values
           -ocom  output complex (pickled Dock.Complex object)
           -autocap    discover and cap chain breaks with NME or ACE residues
           -indi       interior dielectric (default 4)
           -exdi       solvent dielectric (default 80)
           -salt       ionic strength in M (default 0.15)
           -ionrad     ion radius (2.0 A)
           -prbrad     probe radius (1.4 A)
           -bndcon     box boundary condition (see Delphi manual, default 4)
           -scale      grid density in points per Angstroem (default 2.3)
           -perfil     the molecules' largest linear dimension relative to the
                       grid in percent (default 60)

Example Use:

    delphiBinding.py -c input.complex -v -scale 1.2

Or starting from PDB:

    delphiBinding.py -i complex.pdb -r 0 1 -l 2 -autocap -v -scale 2.1
                       
Note: 
     
      -scale and -perfil together determine the size of the grid. The
      recommended value is -scale 2.3. A -perfil of 60% means that, along its
      largest dimension, the solute molecule will still only fill 60% of the
      grid (box). This factor can probably be increased (resulting in a
      smaller grid) for larger complexes.
      
      See Biskit/Dock/DelphiBindingEnergy.py and Biskit/Delphi.py for a more
      detailed discussion of parameters and issues.
      
      The output complex can be opened with Biskit.tools.load(). For example:
      >>> import Biskit.tools as T
      >>> com = T.load( 'delphi.complex' )

      The delphi results can then be accessed from the complex info record 
      'dG_delphi':
      >>> print "free energy of binding: ", com['dG_delphi']['dG_kcal']
      >>> com.model().writePdb('delphi_complex.pdb')
      
      The Biskit script bis.py provides a short cut from the shell:
      ~> bis.py delphi.complex
      ...
      Loaded delphi.complex and put it into variable x.
      >>> x['dG_delphi']['dG_kcal']
      12.34
      
Default options:
"""
    for key, value in options.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)

def inputComplex( options ):
    if 'c' in options:
        return T.load( T.absfile( options['c'] ) )

    m = PDBModel( T.absfile( options['i'] ) )
    
    ## extract rec and lig chains
    rec_chains = T.toIntList( options['r'] )
    lig_chains = T.toIntList( options['l'] )

    rec = m.takeChains( rec_chains )
    lig = m.takeChains( lig_chains )

    ## create Protein complex
    com = Complex( rec, lig )
    return com

def report( com ):
    t = '%4s : %4s %7s %7s %7s %7s %7s' %\
      ('Rec','Lig','dG(kT)','solv','coul','ions','dG(kcal/mol)\n')

    s = '%(rec)4s : %(lig)4s %(dG_kt)7.2f %(solv)7.2f %(coul)7.2f %(ions)7.2f %(dG_kcal)7.2f\n'
    
    d = copy.copy(com['dG_delphi'])
    d.update( {'rec':com.rec_model.pdbCode, 'lig':com.lig_model.pdbCode} )
    
    r = t
    r += s % d
    return r
    
### MAIN ###
############

options = T.cmdDict( {'ocom':'delphi.complex', 'o':'dg_delphi.out' } )

try:

    f_out = T.absfile( options['o'] )
    f_ocom = T.absfile( options['ocom'] )
    
    options['autocap'] = 'autocap' in options
    options['debug'] = 'debug' in options
    options['verbose'] = 'v' in options
    
    if 'nice' in options: options['nice'] = int( options['nice'] )

    for key in ['indi', 'exdi', 'salt', 'ionrad', 'prbrad', 
                 'bndcon', 'scale', 'perfil']:
        if key in options: options[key] = float( options[key] )

    if 'log' in options:
        options['log'] = LogFile( options['log'] )

    ## create a complex
    com = inputComplex( options )

    dg = DelphiBindingEnergy( com, **options )
    r = dg.run()
    
    print "Saving result complex to ", f_ocom
    T.dump( dg.delphicom, f_ocom )
    
    print "Final Result"
    print "============"
    print report( dg.delphicom )
    
    f = open( f_out, 'w' )
    f.write( report( dg.delphicom ) )
    f.close()
    print "energy values written to ", f_out
    
except KeyError, why:
    print 'Insufficient options. Missing: ', (str(why))
    _use( options )

except Exception, why:
    print "There was an error..."
    print T.lastError()
    print T.lastErrorTrace()






    


