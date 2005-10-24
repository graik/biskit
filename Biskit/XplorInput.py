##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## last $Author$
## $Date$

import Biskit.tools as t
import os

class XplorInputError(Exception):
    pass

class XplorInput:
    """
    Create Xplor input file. addFromTemplate() might be usefull for
    non-xplor-stuff, too.
    Note, the file is only flushed to disc when the object is destructed,
    or the flush() method is called!
    """

    def __init__(self, outName, mode='w'):
        """
        outName - Folder to place output files in (String)
        baseName- start of generated filenames, i.e. |1AKZ_|generate.inp
        mode    - open file with this mode, w=override, a=append
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
        """
        try:
            self.fgenerate.write(str + '\n')
        except (IOError):
            t.errWriteln(
                "XPlorInput.append(): Error adding string to xplor input file.")
            t.errWriteln( t.lastError() )


    def _singleParam(self, param, value, indent):
        """
        return single parameter=value line with |indent| leading tabs.
        (Used by blockFromDic)
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
        title   - first line of block
        paramDic- dictionary of type {'param1':value1, 'param2':value2,...}
        indent  - number of tabs to add before each line
        priorityParams  - list of params to write first (in that order) e.g.
                          [param1, param2] to have param1 and param2 come first
        Will result in param1=value1, param2=value2.
        If e.g. value1 is an empty string or None just param is written in that
        line
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

        title    - string, first line
        paramLst - list of strings, i.e. ['nsteps=10', 'eps=1']
        indent   - int, number of tabs indentation

        -> string, block with parameters.
        """
        result = indent*'\t' + title + '\n'  # first line of block
        
        for param in paramLst:
            result = result + (indent+1) * '\t' + param + '\n'
        result = result + indent*'t' + 'end\n'

        return result


    def addBlockFromDic(self, title, paramDic, indent=0):
        """
        Convenience implementation of block. This one directly appends the block
        to the growing file.
        """
        self.add( self.blockFromDic(title, paramDic, indent) )


    def addBlock(self, title, paramLst, indent=0):
        """
        create block and directly add it to growing input file.
        """
        self.add( self.block( title, paramLst, indent) )


    def renameAtom(self, resname, oldName, newName):
        """
        append statement for renaming atoms:
        -> vector do (name=|newName|) (name |oldName| and resname |resname|)
        """
        self.add(
            ("""vector do (name="%(newName)s") (name "%(oldName)s" """+
            """and resname "%(newName)s" )""") % vars() )


    def renameAtoms( self, atomDicLst ):
        """
        Rename several atoms.
        atomDicLst - list of dictionaries ala
                     [{'res':'ALA', 'old':'HT1', 'new':'H1'},{...}]
                     each dictionary describing one renaming task """
        
        for rename in atomDicLst:
            self.renameAtom(rename['res'], rename['old'], rename['new'])


    def renameRes(self, oldName, newName, select=''):
        """
        Append statement for renaming residues:
        oldname - str, current residue name
        newname - str, new one
        select  - optional additional selection for residue, e.g. 'and resid 10'
        -> str, 'vector do (resname=|newName|) (resname |oldName| |select|)' """
        
        self.add(\
            ("""vector do (resname="%(newName)s") (resname "%(oldName)s" """+
             """%(select)s )""") % vars() )
        

    def hbuild(self, selection):
        """
        Append hbuild block. Example for how to use the block() method.
        selection - String with atom selection
        """
        selection_value = "( %s )" % selection
        result = self.blockFromDic('hbuild',{'selection': selection_value,
                                             'print':'on'})
        self.add( result )
        

    def patchSS(self, res1, res2):
        """
        Patch statement for disulfide bond between residue number 1 and 2.
        res1, res2 - {'res':19,'id':'SEGID'}, number and segid for both residues
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
        segid   - segment id
        fpdb    - complete filename of pdb
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
        fTemplate - String, filename for template
        valueDic  - Dictionary, {placeHolder:value}
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
                
            s += "\n  Error:\n  " + t.lastError()

            raise XplorInputError, s 


##################################################
# main function for testing only
##################################################
if __name__ == '__main__':

    test = XplorInput('TEST')
    test.add("\nremarks test generate.inp\n")

    test.addBlockFromDic("minimize powell",{"nsteps":100,"npr":5})
    test.addAmberSegment("id_1", "/home/Bis/super.pdb")
    test.patchSS({'res':23, 'id':'id_1'},{'res':44, 'id':'id_1'} )
    test.hbuild("hydrogen and not known")

    test.renameRes("CYS", "CYX")
    test.renameAtom("ARG","HG1","HG2")
    test.addFromTemplate("/home/Bis/raik/data/tb/xplor/t_template.dat",
                         {'value':'TEST','number':10})
    test.flush()
##    test = 1
