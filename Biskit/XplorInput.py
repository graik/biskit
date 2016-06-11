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
Create Xplor input files.
"""

import Biskit.tools as T
from Biskit import EHandler
import os

class XplorInputError(Exception):
    pass

class XplorInput:
    """
    Create Xplor input file. addFromTemplate() might be usefull for
    non-xplor-stuff, too.

    @note: The file is only flushed to disc when the object is destructed,
    or the flush() method is called!
    """

    def __init__(self, outName, mode='w'):
        """
        @param outName: Folder to place output files in
        @type  outName: str
        @param mode: open file with this mode, w=override, a=append
        @type  mode: str
        """
        self.foutName = os.path.abspath(outName)    # name of new .inp file

        # open for <appending|writing|reading>
        self.fgenerate = open(self.foutName, mode) 


    def __del__(self):
        self.fgenerate.close()


    def flush(self):
        """
        Flush output file (but keep it open).
        """
        self.fgenerate.flush()


    def add(self, str):
        """
        Add String str and line break to xplor input file.

        @param str: string to add to file
        @type  str: str        
        """
        try:
            self.fgenerate.write(str + '\n')
        except (IOError):
            EHandler.error(
                "XPlorInput.append(): Error adding str to xplor input file.")

    def _singleParam(self, param, value, indent):
        """
        return single parameter=value line with |indent| leading tabs.
        (Used by blockFromDic)

        @param param: parameter
        @type  param: str
        @param value: value
        @type  value: str
        @param indent: number of tabs at begining of line
        @type  indent: int

        @return: line vith indent+parameter+value
        @rtype: str
        """
        result = indent*'\t' + param
        try:
            if (value <> None) and (value <> ''):
                result = result + "=" + str(value)
        except:
            pass
        return result + '\n'


    def blockFromDic(self, title, paramDic, indent=0, priorityParams=None):
        """
        Returns block with parameters (as String).
        Use e.g. for minimize or paramter statement.

        @param title: first line of block
        @type  title: str
        @param paramDic: dictionary of type
                         C{ {'param1':value1, 'param2':value2,...} }
        @type  paramDic: dict
        @param indent: number of tabs to add before each line
        @type  indent: int
        @param priorityParams: list of params to write first (in that order)
                               e.g. C{ [param1, param2] } to have param1 and
                               param2 come first
        @type  priorityParams: [str]

        @return: Will result in param1=value1, param2=value2. If e.g. value1
                 is an empty string or None just param is written in that line
        @rtype: str
        """
        result = ""
        if indent > 0:
            # if nested, the first line most likely already has a ta
            result = (indent-1)*'\t'

        result = result + title + '\n'  # first line of block

        # create ordered list of (parameter, value) tuples
        paramItems = []
        # put priority Params first
        if priorityParams <> None:
            for param in priorityParams:
                paramItems = paramItems + [(param, paramDic[ param ])]

        # add all others Params
        for param in paramDic.keys():   # add all others next
            if (param, paramDic[ param ]) not in paramItems:
                paramItems = paramItems + [(param, paramDic[ param ])]
        for paramTuple in paramItems:
            result += self._singleParam(paramTuple[0], paramTuple[1], indent+1)

        result = result + indent*"\t" + "end"   # end of block
        return result


    def block(self, title, paramLst, indent=0):
        """
        As block() but takes parameters as list of strings without resorting
        to the clumpsy dictionary.

        @param title: first line
        @type  title: string
        @param paramLst: list of strings, i.e. C{ ['nsteps=10', 'eps=1'] }
        @type  paramLst: [str]
        @param indent: number of tabs indentation
        @type  indent: int

        @return: block with parameters.
        @rtype: str
        """
        result = indent*'\t' + title + '\n'  # first line of block

        for param in paramLst:
            result = result + (indent+1) * '\t' + param + '\n'
        result = result + indent*'t' + 'end\n'

        return result


    def addBlockFromDic(self, title, paramDic, indent=0):
        """
        Convenience implementation of block. This one directly appends
        the block to the growing file.

        @param title: first line of block
        @type  title: str
        @param paramDic: dictionary of type
                         C{ {'param1':value1, 'param2':value2,...} }
        @type  paramDic: dict
        @param indent: number of tabs to add before each line
        @type  indent: int        
        """
        self.add( self.blockFromDic(title, paramDic, indent) )


    def addBlock(self, title, paramLst, indent=0):
        """
        create block and directly add it to growing input file.

        @param title: first line
        @type  title: string
        @param paramLst: list of strings, i.e. C{ ['nsteps=10', 'eps=1'] }
        @type  paramLst: [str]
        @param indent: number of tabs indentation
        @type  indent: int        
        """
        self.add( self.block( title, paramLst, indent) )


    def renameAtom(self, resname, oldName, newName):
        """
        append statement for renaming atoms:

        @param resname: name of residue to replace atom in  
        @type  resname: str
        @param oldName: atom name to replace
        @type  oldName: str
        @param newName: new atom name
        @type  newName: str

        @return: Line: C{ vector do (name=|newName|) (name |oldName|
                          and resname |resname|) }
        @rtype: str
        """
        self.add(
            ("""vector do (name="%(newName)s") (name "%(oldName)s" """+
             """and resname "%(newName)s" )""") % vars() )


    def renameAtoms( self, atomDicLst ):
        """
        Rename several atoms.

        @param atomDicLst: list of dictionaries wher each dictionary
                           describs one renaming task e.g. ::
                             [{'res':'ALA', 'old':'HT1', 'new':'H1'},{...}]
        @type  atomDicLst: dict
        """ 
        for rename in atomDicLst:
            self.renameAtom(rename['res'], rename['old'], rename['new'])


    def renameRes(self, oldName, newName, select=''):
        """
        Append statement for renaming residues:

        @param oldName: current residue name, to replace
        @type  oldName: str
        @param newName: new residue name
        @type  newName: str
        @param select: optional additional selection for residue,
                       e.g. 'and resid 10'
        @type  select: str

        @return: Line: C{ 'vector do (resname=|newName|)
                           (resname |oldName| |select|)' }
        @rtype: str
        """

        self.add(\
            ("""vector do (resname="%(newName)s") (resname "%(oldName)s" """+
             """%(select)s )""") % vars() )


    def hbuild(self, selection):
        """
        Append hbuild block. Example for how to use the block() method.

        @param selection: String with atom selection
        @type  selection: str
        """
        selection_value = "( %s )" % selection
        result = self.blockFromDic('hbuild',{'selection': selection_value,
                                             'print':'on'})
        self.add( result )


    def patchSS(self, res1, res2):
        """
        Patch statement for disulfide bond between residue number 1 and 2.

        @param res1: dictionary with number and segid for the first residue
                     e.g. C{ {'res':19,'id':'SEGID'} }
        @type  res1: dict
        @param res2: dictionary with number and segid for the second residue
        @type  res2: dict
        """
        select_1 = """( segid "%4s" and resid %s )""" % (res1['id'],res1['res'])
        select_2 = """( segid "%4s" and resid %s )""" % (res2['id'],res2['res'])
        result = self.blockFromDic("patch disu",
                                   {'reference=1': select_1, 'reference=2': select_2} )
        self.add(result)


    def addAmberSegment(self, seg_id, fpdb):
        """
        Return string holding xplor input for segment consisting of one chain.
        Example for nesting blocks with blockFromDic().

        @param seg_id: segment id
        @type  seg_id: str
        @param fpdb: complete filename of pdb
        @type  fpdb: str
        """
        ## read sequence from coordinates
        coord_statement = """coordinates @%(fpdb)s""" % vars()  

        ## create inner block "chain", put coordinates statement at top
        chain_block = self.blockFromDic("chain",\
                                        {'link ppgp head - GLY tail + PRO end':'',
                                         'link ppgg head - GLY tail + GLY end':'',
                                         'link pepp head - * tail + PRO end':'',
                                         'link ppg2 head - GLY tail + * end':'',
                                         'link ppg1 head - * tail + GLY end':'',
                                         'link pept head - * tail + * end':'',
                                         coord_statement:''}, 1, [coord_statement])

        ## create outer block "segment" and nest chain block into it
        result = self.blockFromDic("segment",
                                   {'name':seg_id, chain_block:''}) 

        self.add( result )


    def addFromTemplate(self, fTemplate, valueDic):
        """
        Read template input file with formatstr placeholders, insert values
        from valueDic.

        @param fTemplate: filename for template
        @type  fTemplate: str
        @param valueDic: Dictionary, {placeHolder:value}
        @type  valueDic: dict
        """
        try:

            line = None
            for line in open(fTemplate).readlines():
                line = line.rstrip()               ## get rid of \n

                ## skip template comments starting with '#'
                if len(line) == 0 or line.lstrip()[0] <> '#':
                    self.add( line % valueDic )

        except KeyError, why:
            s =  "Unknown option in template file."
            s += "\n  template file: " + fTemplate
            s += "\n  Template asked for a option called " + str( why[0] )
            s += "\n  template line:\n  " + str(line)
            s += "\n  Please give a value for this option at the command line."
            s += "\n  E.g: pdb2xplor.py -i input.pdb -%s some_value -t..." %\
              str( why[0] )

            raise XplorInputError, s

        except:
            s =  "Error while adding template file."
            s += "\n  template file: " + fTemplate
            s += "\n  template line:\n  " + str(line)
            s += "\n  available arguments:\n"

            for i in valueDic.keys():
                s += "\t%25s\t%s\n" % (i, str( valueDic[i] ) )

            s += "\n  Error:\n  " + T.lastError()

            raise XplorInputError, s 


#############
##  TESTING        
#############

import Biskit.test as BT
import tempfile

class Test( BT.BiskitTest ):
    """
    XPlorer Test
    """

    def prepare(self):
        ## create an temporary input template
        self.f_inp = tempfile.mktemp('_test.inp')
        self.f = open( self.f_inp, 'w')
        self.f.write('Test Template with %(number)i values:\n')
        self.f.write('\n')
        self.f.write('%(value)s\n')
        self.f.close()

        ## temporary output template
        self.f_out_inp = tempfile.mktemp('_test_out.inp')

    def test_XplorInput( self ):
        """XplorInput test"""

        ## write to output
        self.t = XplorInput( self.f_out_inp )
        self.t.add("\nremarks test generate.inp\n")

        self.t.addBlockFromDic("minimize powell",{"nsteps":100,"npr":5})
        self.t.addAmberSegment("id_1", "/home/Bis/super.pdb")
        self.t.patchSS({'res':23, 'id':'id_1'},{'res':44, 'id':'id_1'} )
        self.t.hbuild("hydrogen and not known")

        self.t.renameRes("CYS", "CYX")
        self.t.renameAtom("ARG","HG1","HG2")
        self.t.addFromTemplate( self.f_inp, {'value':'TEST','number':10})
        self.t.flush()

        self.__dict__.update( locals() )

        ## check result
        f = open( self.f_out_inp, 'r')
        self.lines = f.readlines()

        if self.local:
            for l in self.lines:
                print l[:-1]
        f.close()

        self.assertEqual( self.lines[-3], 'Test Template with 10 values:\n' )

    def cleanUp(self):
        T.tryRemove( self.f_inp )
        T.tryRemove( self.f_out_inp )


if __name__ == '__main__':

    BT.localTest()



