from Biskit import Executor,BiskitError
from socketoutparser import parse
from Biskit.DSSP import Dssp
import tempfile
import os

class Socket_Error( BiskitError ):
    pass


class SocketCoil (Executor):
    def __init__(self,model,**kw):
        self.model = model
        
        ## temporary pdb-file
        self.f_pdb = tempfile.mktemp( '_sock.pdb')
        self.f_dssp = tempfile.mktemp( '_sock.dssp')
        self.f_out = tempfile.mktemp( '_sock.out')
        
        print
        print self.f_pdb
        print self.f_dssp
        print self.f_out
        
        Executor.__init__( self, 'socket',
                           args='-f %s -s %s'%(self.f_pdb,self.f_dssp),
                           catch_err=1, **kw )
                           
    
    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.dssp = Dssp(self.model,f_out = self.f_dssp)
        self.dssp.run()
        
        ## Check for output file
        assert (os.path.exists(self.f_dssp)), "Error while creating dssp file."
        
        self.mask  = self.model.maskProtein(standard=True)
        self.model = self.model.compress( self.mask )
        self.model.writePdb( self.f_pdb )


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )
            T.tryRemove( self.f_dssp )
            T.tryRemove( self.f_out )


    def parse_result( self ):
        print "parsing"
        print self.f_out
        return parse(self.f_out)
    
    
    def isFailed(self):
        return False
    
    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        ## Poner en el modelo como profile el registro model.['ccsocket'] ='abcd'
        ## self.model.residues.set( 'ccsocket', profile, mask=self.resMask )
        self.result = self.parse_result()
        
        
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit import PDBModel

class Test(BT.BiskitTest):
    """ Test cases for Coiled Coil Utils"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    def test_execution(self):
        """ testing of execution"""
        p = PDBModel(T.testRoot()+"/coiledcoil/1ZIJ.pdb")
        s = SocketCoil(p)
        s.debug = True
        s.run()
        
if __name__ == '__main__':
    BT.localTest()    
    