##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
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
    hosts: list of hostnames
    """
    for host in hosts:
        try:
            result = addhosts([host])
        except:
            print host,'failed'

def delHosts(hosts):
    return delhosts( hosts )
