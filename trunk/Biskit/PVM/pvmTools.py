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
##
## last $Date$
## last $Author$
"""low-level tools for PVM"""

import pvm
import pypvm_core
import os
import Biskit.hosts as H


def expandLocal( host ):
    """
    Expand name of local(!) host.
    For example, turn 'myhost' into 'myhost.mydomain.com'
    """
    if os.uname()[1].split('.')[0] == host:
        return os.uname()[1]
    
    return host


def addHosts( number=len(H.nodes_all), hosts=H.nodes_all, expand=1 ):
    """
    Add hosts to PVM.
    number   - int, number of requested nodes
               if number == 0, nothing happens
    hosts    - [ str ]
    localise - [1|0], add full address to local node
    """
    if expand:
        hosts = [ expandLocal( h ) for h in hosts ]
    
    if number == None:
        pvm.addHosts( hosts )

    if number > 0:
        pvm.addHosts( hosts[:number+1] )


def countHosts():
    """
    -> int, number of hosts registered with PVM
    """
    return PVM.pypvm_core.config()[0]
