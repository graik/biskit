##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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
## last $Date$
## last $Author$

"""
PVM interface.
"""

## MAX_RECURSION_LIMIT = 100000

from Biskit.PVM.pypvm_core import *

from cPickle import dumps # as _dumps
from cPickle import loads # as _loads


## def dumps(object):
##     import sys

##     limit = sys.getrecursionlimit()

##     while limit < MAX_RECURSION_LIMIT:

##         try:
##             s = _dumps(object)

##             o = _loads(s)

##             if not type(o) == type(object):
##                 print 'dumps: inconsistency; %s, %s' % \
##                       (`type(o)`, `type(object)`)

##             return s

##         except RuntimeError:
##             limit += 1000
##             sys.setrecursionlimit(limit)

##     raise RuntimeError, "Maximum recursion depth exceeded."

## def loads(object):
##     import sys

##     limit = sys.getrecursionlimit()

##     while limit < MAX_RECURSION_LIMIT:

##         try:
##             return _loads(object)

##         except RuntimeError:
##             limit += 1000
##             sys.setrecursionlimit(limit)

##     raise RuntimeError, "Maximum recursion depth exceeded."


def pack(object):
    return pkstr(dumps(object))


def unpack():
    return loads(upkstr())


def pack_and_send(tid, msg_tag, object, encoding = None):
    if encoding is None:
        encoding = data['default']

    return psend_str(encoding, tid, msg_tag, dumps(object))


def addHosts(hosts):
    """
    Add hosts to PVM.
    
    @param hosts: list of hostnames
    @type  hosts: [str]
    """
    for host in hosts:
        try:
            result = addhosts([host])
        except:
            print host,'failed'


def delHosts(hosts):
    """
    Remove hosts from PVM.
    
    @param hosts: list of hostnames
    @type  hosts: [str]
    """
    return delhosts( hosts )
