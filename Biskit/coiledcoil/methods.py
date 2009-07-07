import Biskit.tools as T
import socket_coil as SC
import paircoil as PC
from coiledutils import getHeptad
from coiledcoil import CoiledCoil



"""
Methods overview:

A method is a function which, given a chain sequence, predicts the heptad register
on it ( So for instance socket is not a method of this kind, as it predicts the
heptad using the structure).

Methods are defined in the 'METHODS' structure, and indexed by a method Id.
So it's a dictionary indexed by Id where:

- First entry is the type of tool used to predict. By now it can be:
    . EXE - if it uses an extern binary file (like Paircoil)
    . TABLE - if it uses one table (stored in /data/coiledcoil)
    . NONE - if it's just a dummy method which doesn't predict anything
        (used for manual fulfilling for instance).

- Second entry is the name of the element used to predict. If first is 'EXE' then
  is the name of the binary, if first is 'TABLE' then it's the name of the table 
  file, and so on...

- Third entry is a short description of the method. The shorter the better.

- Fourth is the name of a function which has to be called for rretrieving
  the register. Used for instance to call the wrapper of a binary and 
  post-process its results.

"""

def paircoil_method(method,seq = None):
    """
    Function for retrieving and post processing results from 'Paircoil2'
    
    @param seq: Sequence to extract registers. 
    @type seq: string
    """
    if len(seq) >= 28: ## has to be like the window to work
        pc = PC.Paircoil(seq)
        pc.debug = True
        pc.run()
        pc.result.register
        return method+":"+getHeptad(pc.result.chain,pc.result.register)
    else:
        return ""

"""
Methods description structure.

Special methods:

'Socket', method used by default. Always present.

'Paper' When this method is used, the format must be : 'Paper:register:reference'
    - 'reference' must not contain spaces.

"""
METHODS = {"Pair":("EXE","","Paircoil",paircoil_method),\
           "SPar":("TABLE","SOCKET_antipar_norm","Socket Parallel Score Table"),\
           "SAPar":("TABLE","SOCKET_antipar_norm","Socket Antiparallel Score Table"),\
           "Paper":("NONE","","Paper (reg+reference)"),\
           "Parry":("TABLE","DADParry_scaled","D A Parry score table"),\
           "Default":("NONE","","Not defined method (just for using defaults)")\
           }


priorities = {"Socket":5,## Socket is the reference  \ 
            "Pair":3,\
            "SPar":4,\
            "SAPar":4,\
            "Paper":5,\
            "Parry":3,\
            "Default":2\
            }

"""
Enumeration of methods available
"""
sources = METHODS.keys() 



def getRegisterByMethod( seq = "", method = ""):
    """
    Core function. Given a sequence and a method Id it returns the heptad
    representative for this chain.
    
    @param seq: Sequence for the heptads to be predicted
    @type seq: string
    @param method: Method Id as stored in 'METHODS' structure 
    @type method: string
    
    @return: A string representing the method and the predicted heptad
    @rtype: string
    """
    assert (method in sources),"Not defined method %s."%method
    
    
    if METHODS[method][0] == "TABLE":
        ## The method relies in a score table in the data folder, so
        ## pick this table and calculate
        cc = CoiledCoil(T.dataRoot()+"/coiledcoil/"+METHODS[method][1])
        return method+":"+cc.findHeptads(seq)['best']
        
    elif METHODS[method][0] == "EXE":
        return METHODS[method][3](method,seq)
    
    elif METHODS[method][0] == "NONE":
        return ""
    
    
    return ""
    
def getCoilByMethod(method = ""):
    """
    Creates a CoiledCoil instance suitable for that particular method.
    
    @param method: Method Id as stored in 'METHODS' structure 
    @type method: string
    
    @return: CoiledCoil with correct scoring (table) info for a method.
    @rtype: CoiledCoil
    """
    if method != "" and METHODS[method][0] == "TABLE":
        cc = CoiledCoil(T.dataRoot()+"/coiledcoil/"+METHODS[method][1])
    else:
        cc = CoiledCoil()
        
    return cc 