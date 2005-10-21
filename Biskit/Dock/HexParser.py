## HexParser:
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## Parse output file from hex docking run
##
## last $Author$
## $Date$

import re, sys
import Numeric  ## array
from Biskit.Dock import Complex, ComplexList

from Biskit import PCRModel
import Biskit.tools as t

class HexParser:
    """
    Parse hex result file and extract Complex objects into a dictionary
    indexed by the hex solution number.
    The hex file is only closed when this object is discarded.
    The HexParser is created with a dictionary containing PCRModel objects
    for ligand and receptor indexed by hex model number.
    """

    def __init__(self, hexFile, rec_dic, lig_dic, forceModel=None):
        """
        hexFile - string, name of hex output file
        rec_dic - dictionary, map between model number and PCRModel
                  { 1 : model1, 2 : model2,..} types: { int : PCRModel}
        lig_dic - same for ligand PCRModels
        forceModel - (int, int), force parser to accept model numbers
        """
        self.hex = open( hexFile )
        self.rec_models = rec_dic
        self.lig_models = lig_dic
        self.forceModel = forceModel
        ## pattern for analyzing single line of type "Value: -1.23":
        self.ex_line = re.compile("^(\S+):\s*([-0-9\.\s]*)")
        ## pattern to find one line of transformation matrix
        self.ex_matrix = re.compile("([-0-9]+\.[0-9]+e[-+0-9]+)")


    def __del__(self):
        self.hex.close()


    def nextComplex(self):
        """
        Take list of lines, extract all Hex info about one complex
        (Solution number, hex energy,..) also extract 16 numbers of
        the transformation matrix and put them into 4 by 4 numeric
        array.
        -> Complex, created from that information
        """
        ## get set of lines describing next complex:
        lines = self._nextBlock()
        if lines == None:
            ## EOF
            return None
        ## skip incomplete records
        if len(lines) < 13:
            lines = self._nextBlock()
            
        ## fill info dictionary
        i = {}
        matrix = None
        for l in lines:
            try:
                ## labels has to be in the same order as in hex.out
                m = self.ex_line.search( l )
                if m != None:
                    m = m.groups()
                
                    if m[0] == 'Orientation':
                        i['hex_clst'] = int(m[1])
                    elif m[0] == 'Solution':
                        i['soln'] = int(m[1])
                    elif m[0] == 'ReceptorModel':
                        if self.forceModel:
                            i['model1'] = self.forceModel[0]
                        else:
                            i['model1'] = int(m[1])
                    elif m[0] == 'LigandModel':
                        if self.forceModel:
                            i['model2'] = self.forceModel[1]
                        else:
                            i['model2'] = int(m[1])
                    elif m[0] == 'Bumps':
                        if int(m[1]) != -1:
                            i['bumps'] = int(m[1])
                    elif m[0] == 'ReferenceRMS':
                        i['rms'] = float(m[1])
                    elif m[0] == 'Vshape':
                        if float(m[1]) != 0.0:
                            i['hex_Vshape'] = float(m[1])
                    elif m[0] == 'Vclash':
                        if float(m[1]) != 0.0:
                            i['hex_Vclash'] = float(m[1])
                    elif m[0] == 'Etotal':
                        i['hex_etotal'] = float(m[1])
                    elif m[0] == 'Eshape':
                        i['hex_eshape'] = float(m[1])
                    elif m[0] == 'LigandMatrix':
                        ## get all numbers of matrix as list of strings
                        strings = self.ex_matrix.findall( l )
                        ## convert that to list of floats
                        numbers = []
                        for each in strings:
                            numbers += [float(each)]
                        ## convert that to list of lists of 4 floats each
                        matrix = []  
                        for j in range(0,4):
                            matrix.append( numbers[4*j:4*(j+1)] )
                        ## create 4 by 4 Numeric array from 4 by 4 list
                        matrix = Numeric.array(matrix, 'f')
            except AttributeError:
                print "HexParser.nextComplex(): ",t.lastError()
                            
        ## Create new complex taking PCR models from dictionary
        c = Complex( self.rec_models[ i['model1'] ],
                     self.lig_models[ i['model2'] ],  matrix, i )
        return c
        

    def _nextBlock(self):
        """
        return all lines describing next complex in hex file.
        -> [str, str, str, ...]
        """
        line = self.hex.readline()
        result = []
        while line:
            if line[:4] != "# --":
                ## found matrix line, append it to previous line
                if self.ex_matrix.search(line) != None:
                    result[-1] = result[-1] + line
                else:
                    result += [line]
            else:
                return result
            line = self.hex.readline()

        ## end of file:
        return None


    def parseHex(self):
        """
        Create one Complex Object for each paragraph in hex output
        file.
        -> ComplexList
        """
        complexes = ComplexList()
        c = self.nextComplex()

        while (c <> None):
            complexes.append(c)
            c = self.nextComplex( ) ## look for next cluster

        return complexes

#############################
## TESTING

if __name__ == '__main__':

    
    
    rec_dic = t.Load( t.testRoot() + "/dock/rec/1A2P_model.dic" )
    lig_dic = t.Load( t.testRoot() + "/dock/lig/1A19_model.dic" )

    h = HexParser( t.testRoot() + "/dock/hex/1A2P-1A19_hex.out",
                   rec_dic, lig_dic)
    c_lst = h.parseHex()
    
    print c_lst[1]
