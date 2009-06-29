from Biskit import Executor,BiskitError
from 


class Socket_Error( BiskitError ):
    pass


class SocketCoil (Executor):
    def __init__(self,model,**kw):
        self.model = model
        
        ## temporary pdb-file
        self.f_pdb = tempfile.mktemp( '_sock.pdb')
        self.f_dssp = tempfile.mktemp( '_sock.dssp')
        self.f_out = tempfile.mktemp( '_sock.out')
        
        Executor.__init__( self, 'socketcoil',
                           args='-f %s -s %s'%(self.f_pdb,self.f_dssp),
                           catch_err=1, **kw )
                           
    
    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.dssp = DSSP(self.model,f_out = self.f_dssp)
        self.dssp.run()
        ## Check for output file
        os.path.exists(self.f_dssp)
        self.mask  = self.model.compress( self.model.maskProtein(standard=True))
        self.resmask = self.model.atom2resMask( self.mask )
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


    def parse_result( self ):
        pass
    
    
    def isFailed(self):
        
        return "True if fails"
    
    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        ## Poner en el modelo como profile el registro model.['ccsocket'] ='abcd'
        self.model.residues.set( 'ccsocket', profile, mask=self.resMask )
        self.result = self.parse_result( )
        
        
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
        p = PDBModel(T.testRoot()+"/coiledcoil/pdbs/1VEW.pdb")
        s = SocketCoil(p)
        
if __name__ == '__main__':
    BT.localTest()    
    