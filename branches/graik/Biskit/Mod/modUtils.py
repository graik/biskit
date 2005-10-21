##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$
"""utility funtions for Mod package"""


def parse_tabbed_file( fname ):
    """
    fname - str, file name of key : value mapping
    -> { key : value, }
    """
    f = open( fname )

    result = {}
    for l in f:
        if not l[0] == '#':

            try:
                fname, chain_id = l.split()
                if not len(fname) == 0:
                    result[ fname ] = chain_id
            except:
                fname = l.strip()
                result[ fname ] = ''

    f.close()

    return result
