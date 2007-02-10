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
##
## last $Date$
## last $Author$
"""
low-level tools for PVM
"""

import pvm
import pypvm
import os
import Biskit.hosts as H


def expandLocal( host ):
    """
    Expand name of local(!) host.
    For example, turn 'myhost' into 'myhost.mydomain.com'

    @param host: host name
    @type  host: str

    @return: full host name
    @rtype:  str
    """
    if os.uname()[1].split('.')[0] == host:
        return os.uname()[1]

    return host


def addHosts( number=len(H.nodes_all), hosts=H.nodes_all, expand=1 ):
    """
    Add hosts to PVM.
    
    @param number: number of requested nodes
                   if number == 0, nothing happens
    @type  number: int
    @param hosts: list of host names
    @type  hosts: [str]
    @param expand: add full address to local node
    @type  expand: [1|0]
    """
    if expand:
        hosts = [ expandLocal( h ) for h in hosts ]

    if number == None:
        pvm.addHosts( hosts )

    if number > 0:
        pvm.addHosts( hosts[:number+1] )


def countHosts():
    """
    @return: number of hosts registered with PVM
    @rtype: int
    """
    return PVM.pypvm.config()[0]

## Error codes returned by pvmlib
pvmerrors = {
    0  : 'PvmOk	     Success',
    -2 : 'PvmBadParam     Bad parameter',
    -3 : 'PvmMismatch     Parameter mismatch',
    -4 : 'PvmOverflow     Value too large',
    -5 : 'PvmNoData       End of buffer',
    -6 : 'PvmNoHost       No such host',
    -7 : 'PvmNoFile       No such file',
    -8 : 'PvmDenied       Permission denied',
    -10: 'PvmNoMem        Malloc failed',
    -12: "PvmBadMsg       Can't decode message",
    -14: "PvmSysErr       Can't contact local daemon",
    -15: 'PvmNoBuf        No current buffer',
    -16: 'PvmNoSuchBuf    No such buffer',
    -17: 'PvmNullGroup    Null group name',
    -18: 'PvmDupGroup     Already in group',
    -19: 'PvmNoGroup      No such group',
    -20: 'PvmNotInGroup   Not in group',
    -21: 'PvmNoInst       No such instance',
    -22: 'PvmHostFail     Host failed',
    -23: 'PvmNoParent     No parent task',
    -24: 'PvmNotImpl      Not implemented',
    -25: 'PvmDSysErr      Pvmd system error',
    -26: 'PvmBadVersion   Version mismatch',
    -27: 'PvmOutOfRes     Out of resources',
    -28: 'PvmDupHost      Duplicate host',
    -29: "PvmCantStart    Can't start pvmd",
    -30: 'PvmAlready      Already in progress',
    -31: 'PvmNoTask       No such task',
    -32: 'PvmNotFound     Not Found',
    -33: 'PvmExists       Already exists',
    -34: 'PvmHostrNMstr   Hoster run on non-master host',
    -35: 'PvmParentNotSet Spawning parent set PvmNoSpawnParent',
    -36: "PvmIPLoopback   Master Host's IP is Loopback"
}
    
