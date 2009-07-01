import Biskit.tools as T
import socket_coil as SC
import paircoil as PC



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

def paircoil_method(seq = None):
    """
    Function for retrieving and post processing results from 'Paircoil2'
    
    @param seq: Sequence to extract registers. 
    @type seq: string
    """
    pc = PC.Paircoil(pdb)
    pc.run()
    result =  pc.result
    


"""
Methods description structure.
"""
METHODS = {"Pair":("EXE","","Paircoil",paircoil_method),\
           "SPar":("TABLE","SOCKET_antipar_norm","Socket Parallel Score Table"),\
           "SAPar":("TABLE","SOCKET_antipar_norm","Socket Antiparallel Score Table"),\
           "Paper":("NONE","","Paper (reg+reference)"),\
           "Parry":("TABLE","DADParry_scaled","D A Parry score table"),\
           "Default":("NONE","","Not defined method (just for using defaults)")
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
    """
    assert (not method in sources),"Not defined method."
    
    sequence = pdb.takeChains([chain]).sequence()
    
    if METHODS[method][0] == "TABLE":
        ## The method relies in a score table in the data folder, so
        ## pick this table and calculate
        cc = CoiledCoil(T.dataRoot()+"/"+METHODS[method][1])
        
        return method+":"+cc.findHeptads(seq)['best']
        
    elif METHODS[method][0] == "EXE":
        return method+":"+METHODS[method][3](pdb)
    
    elif METHODS[method][0] == "NONE":
        return sequence[0:7]+":if_paper_put_reference_here"
        ## ideally extract it from a papers library


  