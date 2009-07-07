from Biskit import Executor,BiskitError
from paircoiloutparser import parse
from alignment import PirAlignment
import tempfile
import os


class Paircoil_Error( BiskitError ):
    pass


class Paircoil (Executor):
    def __init__(self,seq,**kw):
        
        self.seq = seq
        
        self.f_fasta = tempfile.mktemp( '_pairc.fasta')
        f_out = tempfile.mktemp( '_pairc.out')
        
        Executor.__init__( self, 'paircoil',
                           args='-win 21 %s %s'%(self.f_fasta,f_out),
                           catch_err=1, f_out = f_out, **kw )
                           
        
        
        #~ print self.f_fasta
        #~ print self.f_out
        
    def prepare( self ):
        """
        Overrides Executor method.
        """
        
        al = PirAlignment()
        al.fasta_b = self.seq
        al.writeFasta(self.f_fasta)


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_fasta )
            T.tryRemove( self.f_out )


    def parse_result( self ):
        
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
        mask  = p.maskProtein(standard=True)
        model = p.compress( mask )
        if self.local:
            print model.sequence()
            
        s = Paircoil(model.sequence())
        s.debug = True
        s.run()
        if self.local:
            print model.sequence()
            print s.result
        
if __name__ == '__main__':
    BT.localTest()    
    