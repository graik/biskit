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

"""
@attention: This class is only used by Prosa.py which should be replaced
by L{Prosa2003}. It will only stay for a while until all transitions are made
to the new version of Prosa.
"""

from UserList import *
import string
from copy import copy

class Table(UserList): 
    def __init__(self, lines = None, lists = None, format = None):
        UserList.__init__(self)
        if lines is not None:
            for line in lines:
                if not format:
                    words = string.split(line)
                else:
                    words = []
                    for item in format:
                        words.append(string.strip(line[item[0]:item[1]]))
                _l = []
                for word in words:
                    try:
                        exec 'item = ' + word
                    except:
                        try:
                            exec 'item = "' + word + '"'
                        except:
                            item = ''
                    if type(item).__name__ in ['builtin_function_or_method',\
                                               'module']:
                        item = word

                    _l.append(item)
                self.data.append(_l)
        else:
            self.data = lists

    def joinColumns(self, left, right, fromRow = -1, toRow = -1, separator = ''):
        if fromRow == -1: fromRow = 0
        if toRow == -1: toRow = len(self.data)
        if left == -1: left = 0
        print fromRow, toRow
        for row in range(fromRow, toRow):
            dummy = copy(right)
            if  dummy == -1: dummy = len(self.data[row])
            item = ''
            for column in range(left, dummy):
                newItem = self.data[row][column]
##                 if not type(newItem).__name__ in ['string','int','float']:
##                     raise 'cannot join types other than string, int, float'
                item = item + separator + str(newItem)
            del self.data[row][left:dummy]
            self.data[row].insert(left,item)

    def __getslice__(self, i,j):
        return Table(lists = self.data[i:j])

    def __getitem__(self, index):
        if type(index).__name__ == 'int':
            return self.data[index]
        else:
            if type(index[0]).__name__ == 'slice':
                lStart = index[0].start
                lStop = index[0].stop
                if lStart is None: lStart = 0
                if lStop is None: lStop = len(self.data)
            else:
                lStart = index[0]
                lStop = lStart + 1

            if type(index[1]).__name__ == 'slice':
                rStart = index[1].start
                rStop = index[1].stop
                if rStart is None: rStart = 0
                if rStop is None: rStop = len(self.data)
            else:
                rStart = index[1]
                rStop = rStart + 1

            columns = self.columns(rStart,rStop)

            return Table(lists = columns[lStart:lStop])

    def __repr__(self):
        return 'table('+UserList.__repr__(self)+')'

    def rows(self, rows):
        try:
            rows[0]
        except:
            rows = (rows,)

        result = []

        for row in rows:
            result.append(self.data[row])

        if len(rows) == 1:
            return asTable(result[0])
        else:
            return asTable(result)

    def columns(self, columns):
        try:
            columns[0]
        except:
            columns = (columns,)

        result = []

        if len(columns) == 1:
            for row in self.data:
                for column in columns:
                    result.append(row[column])
        else:
            for row in self.data:
                col = []
                for column in columns:
                    col.append(row[column])
                result.append(col)

        return asTable(result)

    def deleteRows(self, rows):
        for row in rows:
            del self.data[row]

    def deleteColumns(self, columns):
        for row in range(len( self.data)):
            for column in columns:
                del self.data[row][column]

    def writeToFile(self, name, separator = ' '):
        f = open(name,'w')
        map(lambda row, ff=f, ss=separator: ff.write(string.join(\
            map(lambda item: str(item), row), ss) + '\n'), self.data)
        f.close()


from Scientific.IO import TextFile
def fromFile(fileName, format = None):

    ll = TextFile.TextFile(fileName).readlines()
    return Table(lines = ll, format = format)

def asTable(list):
    return Table(lists = list)



################
## empty test ##
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test"""

    TAGS = [BT.EXE, BT.OLD]

