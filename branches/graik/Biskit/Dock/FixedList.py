#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

class FixedListError( Exception ):
    pass

class FixedList( list ):
    """
    List that blocks any modifications. Implementing lists must override
    fixed() to allow adding of items under some circumstances (like during
    creation of the list).
    """

    def fixed(self):
        """-> true, if list is fixed and modifications are prohibited"""
        return 1

    def __stop(self):
        if self.fixed():
            raise GuardedListError("Attempt to modify fixed list.")
    
    def __setitem__(self, i, v ):
        self.__stop()
        super( GuardedList, self).__setitem__( i, v )

    def __setslice__(self, i, j, v ):
        self.__stop()
        super( GuardedList, self).__setslice__( i, j, v )

    def __add__( self, lst ):
        self.__stop()
        super( GuardedList, self).__add__( lst )

    def __iadd__( self, lst ):
        self.__stop()
        super( GuardedList, self).__iadd__( lst )        

    def extend( self, lst ):
        self.__stop()
        super( GuardedList, self).extend( lst )        

    def append( self, v ):
        self.__stop()
        super( GuardedList, self).append( v )        

    def __delitem__(self, i ):
        self.__stop()
        super( GuardedList, self).__delitem__( i )
