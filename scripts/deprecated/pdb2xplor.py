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


from Scientific.Geometry import Vector

from Biskit.tools import *

from Biskit import LogFile
from Biskit import XplorInput

from Biskit.deprecated.ChainCleaner import ChainCleaner
from Biskit.deprecated.ChainSeparator import ChainSeparator
from Biskit.deprecated.ChainWriter import ChainWriter

from Biskit.Errors import XplorInputError

from Biskit import Executor
from Biskit import settings
from Biskit import Pymoler

## from numpy.oldnumeric import *
## migration not tested (!!):
from Biskit.oldnumeric import *

from string import *
import sys        # sys.exc_type, os.abspath
import re, os, os.path
import commands   # getstatusoutput()


###############################
## create XPlor input script ##

class Xplor (XplorInput):

    def __init__(self, path, cleaner, fheader, fsegment, ftail,
                 amber=0, extras={}):
        """
        path - string, output path for generate.inp
        cleaner - ChainCleaner
        fheader  - string, template for head of generate.inp (file name)
        fsegment - string, template for segment statement (file name)
        ftail    - string, template for end of generate.inp (file name)
        amber    - boolean, switch on CYS -> CYX renaming and add amber patches
        """
        XplorInput.__init__(self, absfile(path) + '/'
                            + cleaner.pdbname + "_generate.inp" )

        self.cleaner = cleaner
        self.log = cleaner.log

        self.path = absfile( path )
        self.path += '/'

        self.header_template = absfile( fheader )
        self.segment_template = absfile( fsegment )
        self.tail_template = absfile( ftail )
        self.amber = amber

        ## announce variables that will be known when creating each segment
        self.amber_patch = ''
        self.segment_id = ''
        self.segment_pdb = ''

        ## pre-defined stuff
        self.ss_cutoff = 4           # threshold of S-S distance for S-S bond
        self.data_root = dataRoot() # folder containing xplor/templates
        self.pdbcode = cleaner.pdbname[:4]
        self.outname = absfile( self.path + self.pdbcode )

        ## add extra fields or override existing ones
        self.__dict__.update( extras )

        self.log.add("\nPreparing XPlor generate.inp file:\n" + 30 * '-' +'\n')

        ## read all segments from cleaner
        self.chains = []
        chain = cleaner.next()
        while chain <> None:
            self.chains = self.chains + [chain]
            chain = cleaner.next()

        ## rename CYS -> CYX for amber FF, only
        if amber:
            self.renameCyx()

        self.log.add("Writing pdb file for each segment.")
        self.writeChains()


    def _getSS(self, res=['CYS','CYX']):
        """
        Get all S-S bonds as list of tuples (each having two residues)
        res - set of string
        """

        ## Create list with all occurances of atomType in residueType
        residueList = []
        for chain in self.chains:

            # keep track of residue's position in residue_list of chain
            resIndex = 0
            for residue in chain.residues:

                try:
                    if residue.name in res: #'CYS':
                        vectorCoord = Vector(residue['SG'].position.array)
                        residueList = residueList +\
                                    [{'res':residue.number, 'id':chain.segment_id,
                                      'xyz':vectorCoord, 'index':resIndex},]
                except:
                    self.log.add(
                        "_getSS(): Error while looking for S-S bonds: ")
                    self.log.add("residue: " + str(residue) )
                    self.log.add("Error: " + lastError() + '\n' )
                    errWriteln("Error in _getSS. See log for details")

                resIndex = resIndex + 1

        ## Calculate distances between all entries in that list
        contactList = []
        counter1 = 0
        while counter1 < len(residueList):
            coord1=residueList[counter1]['xyz']
            counter2 = counter1+1
            while counter2 < len(residueList):
                coord2=residueList[counter2]['xyz']
                if (coord1 - coord2).length() < self.ss_cutoff:
                    contactList = contactList +\
                                [(residueList[counter1], residueList[counter2]),]
                counter2 =counter2 +1
            counter1 = counter1 + 1

        return contactList


    def renameCyx(self):
        """
        Rename all Cys involved in S-S bonds to CYX
        """
        contactList = self._getSS()

        for ssTuple in contactList:
            chainID1 = ssTuple[0]['id']
            chainID2 = ssTuple[1]['id']

            # residue's position in residue_list
            res1, res2 = (ssTuple[0]['index'], ssTuple[1]['index']) 

            for chain in self.chains:
                if chain.segment_id == chainID1:
                    chain.residues[res1].name = 'CYX'
                    self.log.add(
                        "\nRenamed residue %i of segment %s from CYS to CYX."
                        % (res1, chainID1))

                if chain.segment_id == chainID2:
                    chain.residues[res2].name = 'CYX'
                    self.log.add(
                        "\nRenamed residue %i of segment %s from CYS to CYX."
                        % (res2, chainID2))


    def _patchAllSS(self):
        """
        Add patch for each S-S bond. Assumes 'CYS' as residue name!
        """
        # prepare patch statements for SS-bonds
        contactList = self._getSS()
        for ssTuple in contactList:

            try:
                self.patchSS(ssTuple[0], ssTuple[1])
                self.log.add(
                    "\nPatched S-S bond between "+str(ssTuple[0]['res'])+"/"\
                    +ssTuple[0]['id']+ " and "+str(ssTuple[1]['res'])+"/"\
                    +ssTuple[1]['id'])
            except:
                self.log.add( "Error while patching S-S bond." )
                self.log.add( "residues: " + str(ssTuple) )
                self.log.add( "Error: " + lastError() + '\n' )
                errWriteln( "Error while patching S-S. See log for details.")


    def _getTerPatch(self, patchName, ref, resnumber, chain):
        """
        Return patch for terminal residue as string.
        """
        try:
            return """patch %4s refe=%s=( segid "%4s" and resid "%s" ) end\n""" \
                   % (patchName, ref, chain.segment_id, str(resnumber))
        except:
            self.log.add("Error while preparing terminal patch: " + patchName + ":")
            self.log.add( lastError() + "\n")
            errWriteln( "Patching error. See log for details. ")


    def _patchTermini(self, chain):
        """
        Return all amber patches for both termini as string
        """
        ## prepare patch for N- and C-terminus
        nTerResNum = str(chain.residues[0].number)
        n2TerResNum = str(chain.residues[1].number)
        cTerResNum = str(chain.residues[-1].number)
        nTerPatch = 'N' + chain.residues[0].name  # First residue specific N-terminus patch
        cTerPatch = 'C' + chain.residues[-1].name  # C-terminus residue specific patch

        # Second residue specific N-terminus patch
        if nTerPatch == 'NPRO':
            nTerPatch2 = 'PROP'
        if nTerPatch == 'NGLY':
            nTerPatch2 = 'GLYP'
        else:
            nTerPatch2 = 'NTER'

        patch = self._getTerPatch(nTerPatch, "nil", nTerResNum, chain)
        patch = patch + self._getTerPatch(nTerPatch2, '"+"', nTerResNum, chain)
        patch = patch + """patch NTR2 refe="-"=( segid "%4s" and resid "%s" ) refe="+"=( segid "%4s" and resid "%s") end\n\n""" %\
              (chain.segment_id, nTerResNum, chain.segment_id, n2TerResNum)

        patch = patch + self._getTerPatch(cTerPatch, "nil", cTerResNum, chain)
        patch = patch + self._getTerPatch("CTER", '"-"', cTerResNum, chain)

        return patch


    def writeChains(self):
        writer = ChainWriter(self.path )
        for chain in self.chains:
            writer.writeChain(chain)


    def generateInp(self):

        try:
            self.log.add("\nAssembling input file from template files.")
            self.log.add("-"*40 )
            self.log.add("  Header template : " + self.header_template)
            self.log.add("  Segment template: " + self.segment_template)
            self.log.add("  Tail template   : " + self.tail_template)
            self.log.add("\nThese options (place holders) are available:")

            for k in self.__dict__.keys():
                self.log.add("\t%25s\t%s" % (k, str(self.__dict__[k])[:65] ))

            self.addFromTemplate( self.header_template, self.__dict__ )

            for chain in self.chains:
                if self.amber:
                    self.amber_patch = self._patchTermini( chain )
                self.segment_id = chain.segment_id
                self.segment_pdb = absfile( self.path + chain.segment_id +\
                                            "_seg.PDB")

                self.addFromTemplate( self.segment_template, self.__dict__ )

            self._patchAllSS()

            self.addFromTemplate( self.tail_template, self.__dict__ )

            self.flush()

        except XplorInputError, msg:
            errWriteln( "Error while generating the input file. "+\
                        "See log for details.")
            self.log.add( str( msg ) )


##################################################
# Usage and default parameters
##################################################

def _defaultOptions():
    return {
        'o':'.',
        't':dataRoot() + '/xplor/templates/param19',
        'c':'0',
        'thickness':'9',
        'pw_dist':'3.5',
        'water_diam':'2.4'}

def _use(options):
    print """
pdb2xplor:  Create xplor generate.inp from PDB. Chains are separated,
            S-S bonds patched, missing atoms added. Optionally, chain
            breaks are capped with ACE and NME residues. To cap also the
            start or end of chains, you need to add the CA atom of a ACE
            (N') or NME (C') residue with 0.0 0.0 0.0 coordinates to the
            chain ends in the input PDB (using the same segid as the rest
            of the chain).

            NOTE!! The pdb file name has to be 5 characters long (in
                   addition to the .pdb extension) and start with a
                   number.
Syntax:
pdb2xplor -i |pdb_input| [-o |output_path| -c |chain_id_offset| -a -cap
          -t |template_folder| -h |header_template| -s |segment_template|
          -e |end_template| -thickness |solvation_layer|
          -cmask |mask| -exe -view -x |file_with_extra_options|]
or 
pdb2xplor -x |file_with_options|

          -a        .. prepare amber FF (CYS->CYX renaming, terminal patches)
          -h,s,e    .. the three template files
          -t        .. folder containing head.inp, segment.inp, end.inp which
                       will be used as templates if -h,-s, or -e are not given
          -c        .. start chain labeling at position c in the alphabet.
                       (i.e. -c 3 means, chains are labeled D, E, etc.)
          -cap      .. add ACE and NME to N- and C-terminal of chain breaks
          -thickness ..the thickness of the solvation layer in Angstrom
          -cmask    .. chain mask for overriding the default sequence identity
                       based cleaning (e.g. 1 0 0 1 0 0 )
          -exe      .. also execute XPLOR, write log to file
          -view     .. show the cleaned pdb file in PyMol

Default values:
    """
    for key in options.keys():
        print '\t', key, '\t', options[key]
    print "pdb2xplor --help for additional info."
    print

    sys.exit()

def _help():
    print """
The xplor input file will be assembled from 3 template files. The template
files should be independent of the particular PDB and should instead contain
place holders which pdb2xplor will replace by actual file names, numbers, etc.
Place holders look like that:
     %(segment_pdb)s  .. means, insert value of variable segment_id as string
All the variables of the Xplor class (see Xplor.__init__()) can be adressed
this way. The available variables are listed in the log file. Some variables,
like segment_pdb, amber_patch, segment_id will only have meaningfull values in
a segment template.

pdb2xplor combines the templates as follows:

one header_template
+ (one segment template per segment)
+ disulfide patches (generated without template)
+ one tail template

the most relevant variables are:

for header:   data_root    .. Biskit/data folder (xplor top and param files)
for segment:  segment_id   .. segid of currently processed segment
              segment_pdb  .. file name of segment pdb (generated)
              amber_patch  .. terminal patches for amber ff (generated)
for tail:     pdbcode      .. first 4 letters of input pdb file name
              outname      .. suggested file name for output pdb and psf
                              (with absolute path but w/o '.pdb' or '.psf')
              path         .. output path (specified with option -o)

For hackers:
All command line options are also available as variables (e.g. i, o, t).
Even more, you can invent any command line option (e.g. -temperature 298)
which will then be available as variable. Taken the example you could
add a place holder %(temperature)i to your template file which would be
translated to 298 (or whatever you specify).

For hackers++:
With option -x you can specify a file containing variable - value pairs ala:
temperature    298   # the temperature in K
steps 100            !! minimization steps

Give one pair per line, only the first 2 words are parsed.

"""
    sys.exit()

##################################################
# main function
##################################################

def main(options):

##     ## get extra options from external file
##     if options.has_key('x'):
##         options.update( _parseExternalOptions( options['x'] ) )

    fname = options['i']            # input pdb file
    outPath = options['o']

    chainMask = options.get('cmask', None)
    if chainMask: chainMask = toIntList(chainMask)

    try:
        if options.has_key('h'):
            fheader = options['h']
        else:
            fheader = options['t'] + '/head.inp'
        if options.has_key('s'):
            fsegment = options['s']
        else:
            fsegment = options['t'] + '/segment.inp'
        if options.has_key('e'):
            ftail = options['e']
        else:
            ftail = options['t'] + '/end.inp'
    except:
        errWriteln("You have to specify either all three template files or\n"+
                   "(option -t) a folder that contains head.inp, segment.inp, end.inp.\n" +
                   "If both -t and -h, -s, or -e are present, files specified with\n" +
                   "-h, -s, -e are preferred.")

    ## check that input pdb file name (stripped of its path and extension)
    ## is at least 5 characters long and starts with a number
    ## Prompt user for renaming of the file.
    name = stripFilename(fname)

    ## the filename is OK, but it has the same name as used for the
    ## XPLOR output therefore we append _original to the input file
    if toInt( name[0] )== None or len(name)==4 :
        print "##### NOTE: ######"
        print "The pdb file name you gave is OK, but exactly"
        print "the same name will be used for the output from"
        print "the XPLOR script. Therefore '_original' will be"
        print "appended to the input pdb file you gave."
        base, ext = absfile(fname).split('.')
        new_fname = base +'_original.' + ext

        os.rename( fname, new_fname )
        fname = new_fname

    if toInt( name[0] )== None or len(name)<4 :
        print "##### WARNING: ######"
        print "The pdb file name you gave is either shorter "
        print "than 4 characters or it doesn't start with a number."
        print "This will cause the X-PLOR job to fail"
        msg = "Do you want to rename the file (Y/N)?"

        if upper( raw_input( msg ) ) == 'Y':
            new_name = upper( raw_input( "Give a new filname:" ) )

            ## check that the new name is OK
            if toInt( new_name[0] )== None or len(new_name)<4 :
                raise StandardError, \
                      'You gave an incorrect new filename. Exiting'

            ## create link
            else:
                new_file = '%s/%s.%s'%(os.path.dirname(absfile(fname)),
                                       new_name,
                                       fname.split('.')[-1])
                os.link(fname, new_file)
                print 'Link from %s to %s created sucessfully'\
                      %(absfile(fname), new_file)
                fname = new_file

    ## switch on Amber specialities ?
    amber = options.has_key('a')

    ## cap N- and C-term of chain breaks?
    capBreaks = options.has_key('cap')

    cleaner = ChainCleaner(
        ChainSeparator(fname, outPath,
                       int(options['c']),
                       capBreaks=capBreaks,
                       chainMask=chainMask) )

    # initialize with output path and base file name for generate.inp file
    xplorer = Xplor(outPath, cleaner, fheader, fsegment, ftail, amber, extras=options )

    xplorer.generateInp()

    ## run X-Plor
    if options.has_key('exe'):

        out, error, returncode = Executor( 'xplor' , strict=0,
                                           f_in=xplorer.cleaner.pdbname + "_generate.inp",
                                           f_out=xplorer.cleaner.pdbname + '_generate.log').run()
    ## Show structure in pymol
    if options.has_key('view'):

        pm = Pymoler( )
        mname = pm.addPdb( xplorer.outname + '.pdb' )
        pm.add('select xray-wat, segi 1XWW')
        pm.add('select added-wat, segi 1WWW')
        pm.add('select hydrogens, elem H')
        pm.add('hide everything, xray-wat OR added-wat OR hydrogens')
        pm.add('select none')

        colors = [ hex2rgb(c, str) for c in hexColors( len(xplorer.chains) )]

        i=0
        for c in xplorer.chains:
            print colors[i], c.segment_id
            pm.add( 'set_color col_%i, %s'%( i, colors[i] ) )
            pm.add('color col_%i, segi %s and elem c'%( i, c.segment_id ) )
            i += 1

        pm.add('zoom all')
        pm.show()



if __name__ == '__main__':

    options = cmdDict( _defaultOptions() )

    if len(sys.argv) < 2:
        _use(options)
    if sys.argv[1] == '-?' or sys.argv[1] == '--help':
        _help()
    else:
        main(options)
