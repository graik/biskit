
import Biskit.tools as T



"""
Methods overview:



"""
def paircoil_method(pdb = None):
    pass    
    
def socket_method(pdb = None):
    pass  


METHODS = {"Pair":("EXE","","Paircoil",paircoil_method),\
           "SPar":("TABLE","SOCKET_antipar_norm","Socket Parallel Score Table"),\
           "SAPar":("TABLE","SOCKET_antipar_norm","Socket Antiparallel Score Table"),\
           "Paper":("NONE","","Paper (reg+reference)"),\
           "Parry":("TABLE","DADParry_scaled","D A Parry score table"),\
           "Sock":("EXE","","Socket",socket_method),\
           "Default":("NONE","","Not defined method (just for using defaults)")
           }



sources = METHODS.keys()



def getRegisterByMethod( pdb = None, chain = 0, method = ""):
    assert (not method in METHODS),"Not defined method."
    
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


  