#!/usr/bin/env python3
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
## Prepare ensemble MD with amber
## requires template_pme_ensemble folder

import os.path
import os
import numpy as N

from biskit.tools import *
from biskit import PDBModel

host_list = ['localhost'] * 11

def defaultOptions():
    return {'template':dataRoot() + '/amber/template_pme_ensemble_amber9',
            'parm':'',
            'crd':'',
            'out':'.',
            'nb_nodes':1,
            'nodes_eq':host_list[0],
            'nodes_prod':host_list[1:],
            'n_members':10,
            'dt':0.002,
            'n_steps':500000,
            'ntwx':500,
            'ntwv':500}

def syntax( options ):
    print("""
    Syntax:
    amber_ensembleMD.py -parm |parm_file| -crd |crd_file| -out |result_folder|
                        -pdb |0_pdb_file|
                        [ -nb_nodes   |n_nodes_per_host|
                          -template   |template_folder|
                          -nodes_eq   |2 hosts for minimiz.|
                          -nodes_prod |10 hosts for prod.|
                          -n_members  |10|
                          -rseed      |int_random_seed|
                          -dt         |production_time_step|
                          -n_steps    |production_step_number|
                          -ntwx       |production_coordinate_writing_interval|
                          -ntwv       |production_velocities_writing_interval|

                          -place_holder1 |value| -place_holder2 |value| ..]

    The script prepares a new folder |out| with all the input and start files
    to run multiple copies of an Amber MD. All input/start files/folders are
    copied from a template folder. Template files and folders ending in 'xx'
    are recreated |n_members| times. Strings ala '%(place_holder)s' in any
    template file in any template folder are replaced by the value of
    self.|place_holder| which can be given at the command line. If
    self.place_holder contains a list of values, each item is only used once
    (e.g. nodes_prod or nodes_eq ).
       
    Requirements: -$AMBERHOME must be set
                  -LAM environment must be set up in .cshrc or .zshenv or etc.
                  -start_eq must be run from first host in nodes_eq.dat !

    Default options:""")
    for key in options.keys():
        print("\t-",key, "\t",options[key])

    sys.exit(0)


class AmberError( Exception ):
    pass

class AmberController:

    def __init__( self, **options ):
        """
        options must contain:
        template - folder with all template folders and files
        parm     - file name (path) to top.parm
        crd      - file name (path) to starting crd
        out      - result output folder
        n_members- number of ensemble copies
        nb_nodes - number of nodes to use per host
        
        options CAN contain:
        lenres   - number of residues to consider as solute [all non-water]
        rseed    - int, number to add to random seed
        
        dt       - float, time step for production run [0.002]
        n_steps  - int, number of MD steps for production run (nstlim) [500000]
        ntwx     - int, write interval for coordinates [500]
        ntwv     - int, write interval for velocities  [500]
        """
        ## take over parameters into own name space and clean them
        self.__dict__.update( options )
        
        self.template = absfile( self.template, resolveLinks=0 )
        self.out = absfile( self.out, resolveLinks=0 )
        self.n_members = int( self.n_members )
        self.nb_nodes = int( self.nb_nodes )
        self.lenres = int( options.get( 'lenres', 0 ) )
        self.lenatoms = 0

        ## create pool of random seeds, add int if requested
        self.rand_seed = list(range(0, self.n_members*200000, 100))
        if getattr( self, 'rseed', None) is not None:
            self.rand_seed = [ i + int(self.rseed) for i in self.rand_seed ]

        ## get some default MD options 
        self.dt      = float( options.get('dt', '0.002' ) )
        self.n_steps = int( options.get('n_steps', '500000' ) )
        self.ntwx    = int( options.get('ntwx', '500') )
        self.ntwv    = int( options.get('ntwv', '500') )


    def deleteVersionControl(self, path):
        """
        Remove all version control folders in path
        @param path: str, parent folder
        """
        hits = ['CVS', '.svn', '.git']

        for root, dirs, files in os.walk(path, topdown=True):
            for name in files:
                if name in hits:
                    os.remove(os.path.join(root, name))
            for name in dirs:
                if name in hits:
                    tryRemove(os.path.join(root, name), verbose=True, tree=True)


    def prepareFolder(self):

        for f in (self.out, self.template, self.parm, self.crd):
            if not os.path.exists( f ):
                raise AmberError('%s not found.' % f )

        cmd = 'cp -d -r --no-preserve=ownership %s/* %s/' % (self.template, self.out)

        flushPrint('copying template files and folders...')
        os.system( cmd )
        flushPrint('done.\n')

        flushPrint('removing svn / CVS / git files...')
        self.deleteVersionControl(self.out)
        flushPrint('done.\n')

        ## ensure that links are relative to out folder
        if self.out != os.getcwd():
            self.parm = absfile( self.parm )
            self.crd  = absfile( self.crd )
            self.pdb  = absfile( self.pdb )

        os.symlink( self.parm, self.out+'/top.parm' )
        os.symlink( self.crd,  self.out+'/0.crd' )
        os.symlink( self.pdb,  self.out+'/0.pdb' )

        ## number of protein residues into self.reslen
        self.parseReference( self.pdb, dry_out= self.out+'/0_dry.pdb' )

        ## create ensemble folders/ files
        flushPrint('creating ensemble copies of all files/folders with "xx"..')
        files = os.listdir( self.out )
        for f in files:
            result = self.createEnsembleFile( self.out + '/' + f )

            if result:
                os.system( 'rm -rf ' + self.out + '/' + f)

        flushPrint('done.\n')

        ## fill in place holders from own __dict__
        flushPrint('replace place holders in all files...')
        self.completeFile( self.out )
        flushPrint('done.\n')

        ## run prepare scripts in each folder (creating links)
        flushPrint('running PREPARE.* in each folder...')
        self.runPREPARE( self.out )
        flushPrint('done.\n')
        

    def parseReference(self, fpdb, dry_out=None ):
        flushPrint("parsing "+fpdb+"...")
        m = PDBModel( fpdb )
        
        solute_res = m.atom2resMask( N.logical_not( m.maskSolvent() )  )
        
        self.lenres = self.lenres or N.sum( solute_res)
        assert isinstance( self.lenres, N.integer )
        
        self.lenatoms = len( m ) - N.sum( m.maskH2O() )
        assert isinstance( self.lenatoms, N.integer)

        if dry_out:
            m.remove( m.maskH2O() )
            m.writePdb( dry_out )
        flushPrint('done.\n')
        

    def createEnsembleFile(self, fname ):
        """create n_members copies of a file or folder"""

        if fname.find('xx') != -1:

            for i in range(0, self.n_members):
                f = fname.replace('xx','%02i' % i )

                cmd = 'cp -d -r %s %s' % (fname, f)

                os.system(cmd)

            return 1
        
        return 0

    ## iterate over all files in a folder,
    ## replace any place holder by entries of self.__dict__,
    ## or, if the __dict__ contains a list for a given place holder,
    ## insert a unique entry of that list (and remove it from the list)

    def runPREPARE(self, fname):
        """
        recursively run all PREPARE.csh scripts in each folder and
        delete them afterwards.
        """
        if os.path.isdir( fname ):
            files = os.listdir( fname )

            old_wd = os.getcwd()
            os.chdir( fname )

            for f in files:
                self.runPREPARE( absfile( fname+'/'+f, 0) )

            os.chdir( old_wd )
                
        else:
            if stripFilename(fname) == 'PREPARE':
                try:
                    r = os.system( fname )
                    if r == 0:
                        os.remove( fname )
                except:
                    EHandler.error("Couldn't execute " + absfile( fname, 0 ) )
        

    def completeFile(self, fname ):
        """rekursively enter directories and fill out place holders in
        each file."""
        if os.path.isdir( fname ):

            if stripFilename( fname ) not in ['CVS','.svn']:

                files = os.listdir( fname )
                for f in files:
                    self.completeFile( absfile( fname+'/'+f, resolveLinks=0) )
                
        else:
            if not os.path.islink( fname ):
                if fname[-9:] == '_template':
                    os.rename( fname, fname[:-9] )
                self.fillTemplate( fname, fname, self.__dict__ )
            

    def __formatWithLists(self, str, dic ):

        ndic = {}
        for key, item in list(dic.items()):
            if type( item ) is list and str.find('%('+key+')') != -1:
                ndic[ key ] = item[0]
                del item[0]
            else:
                ndic[ key ] = item

        return str % ndic


    def fillTemplate(self, fTemplate, fout, valueDic):
        """
        Read template input file with formatstr placeholders, insert values
        from valueDic. Write it back.
        fTemplate - String, filename for template
        valueDic  - Dictionary, {placeHolder:value}
        """
        try:

            result = []
            line = None
            for line in open( fTemplate ):
                result += self.__formatWithLists( line, valueDic )

            f = open( fout, 'w' )
            f.writelines( result )
            f.close()

        except KeyError as why:
            s =  "Unknown option in template file."
            s += "\n  template file: " + fTemplate
            s += "\n  Template asked for a option called " + str( why[0] )
            s += "\n  template line:\n  " + repr(line)
            s += "\n  Please give a value for this option at the command line."
            s += "\n  E.g: runAmber.py -parm ... -%s some_value -out..." %\
                 str( why[0] )

            raise AmberError(s)
        
        except:
            s =  "Error while adding template file."
            s += "\n  template file: " + fTemplate
            s += "\n  output file: " + fout
            s += "\n  template line:\n  " + str(line)
            s += "\n  available arguments:\n"

            for i in valueDic.keys():
                v = str( valueDic[i] )
                if len( v ) > 120: v = v[:120]
                s += "\t%25s\t%s\n" % (i, v)
                
            s += "\n  Error:\n  " + lastError()

            raise AmberError(s) 

##########
## MAIN ##
##########

if __name__ == '__main__':

    if len( sys.argv ) < 2:
        syntax( defaultOptions() )

    options = cmdDict( defaultOptions() )

    a = AmberController( **options )

    a.prepareFolder()
