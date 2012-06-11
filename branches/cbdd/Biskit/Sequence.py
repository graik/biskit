import Biskit as B
import os
import sys

class SequenceAlignment():
    """

    """
    def __init__(self,referenceModel, pdbList, alignment,seqDict,submask,proteins):
        self.reference=referenceModel
        self.refMol=B.PDBModel(self.reference)
        self.pdbList=pdbList
        self.alignment=alignment
        self.seqDict=seqDict
        self.submask=submask
        self.proteins=proteins

    def superposeAll(self,outputPath="."+os.sep):

        self.refMask=self.refMol.maskFrom("residue_number",self.submask)*self.refMol.maskCA()
        self.resIdxs=[npy.where(self.alignment[0,:]==selResN) for selResN in self.submask]

        i=1
        for c,pdbCode in enumerate(self.pdbList):
            if(pdbCode!=self.reference):
                curMol=self.proteins[i-1]
                submaskCur=[int(self.alignment[i,curel]) for curel in self.resIdxs]
                curMask=curMol.maskFrom("residue_number",submaskCur)*curMol.maskCA()
                if npy.sum(curMask) == npy.sum(self.refMask):
                    (newPos,rmsd)=geom.superpose3D(curMol.xyz, self.refMol.xyz, weights=None,refmask=curMask,targetmask=self.refMask)
                    curMol.xyz=newPos
                    curMol.writePdb(outputPath+"aligned_"+pdbCode)
                i+=1

def alignAndSuperposePdbList(reference,pdbList, submask):
    """
    Aligns and superposes a list of PDB files to a given reference structure using
    the Kabsch structure based alignment algorithm.
    A submask can be specified if the alignment has to be done on a special part
    of the structure. Sequence alignment is performed using muscle.
    TODO --> clean up this and embed muscle better

    @parameter 1 :  reference - PDBModel to take as reference for the structural superimposition
    @type 1 :       PDBModel (or string -> filename)
    ---
    @parameter 2 :  pdbList - Python list of paths to PDB files
    @type 2 :       list
    ---
    @parameter 3 :  submask - A mask indicating
    """

    #submask=npy.array(submask)
    #print submask

    #get sequences for each PDB File
    seqList=[]
    out=open("tmp.fas","w")
    i=0
    print "parsing structures"
    inputSeq=""
    refMol=B.PDBModel(reference)

    #print npy.sum(refMask)
    out.write(">"+reference+"\n")
    inputSeq+=">"+reference+"\n"
    out.write((refMol.compress(refMol.maskProtein())).sequence()+"\n")
    inputSeq+=(refMol.compress(refMol.maskProtein())).sequence()+"\n"

    for pdbPath in pdbList:
        #print "%s - done %f %%"%(pdbPath,float(i)/len(pdbList)*100.)
        curMol=B.PDBModel(pdbPath)
        #seqList.append((curMol.compress(curMol.maskProtein())).sequence())
        out.write(">"+pdbPath+"\n")
        inputSeq+=">"+pdbPath+"\n"
        out.write((curMol.compress(curMol.maskProtein())).sequence()+"\n")
        inputSeq+=(curMol.compress(curMol.maskProtein())).sequence()+"\n"
        i+=1
    out.close()
    print "performing sequence alignment"
    os.system("muscle -quiet -in tmp.fas -out alignment.ala")
    o=open("alignment.ala","r")
    alignment=o.readlines()
    o.close()
    refseq=""
    seqDict={}
    for line in alignment:
        if line[0]==">":
            curCode=line[1:].strip()
            seqDict[curCode]=""
        else :
            seqDict[curCode]+=line.strip()

    seqAlignment=npy.zeros([len(pdbList)+1,len(seqDict[reference])],dtype="int32")
    rn= npy.unique(refMol.compress(refMol.maskProtein())["residue_number"])
    i=0
    for c,el in enumerate(seqDict[reference]):
        if el!="-":
            seqAlignment[0,c]=rn[i]
            i+=1
        else :
            seqAlignment[0,c]=-1

    #print seqAlignment[0,:]
    j=1
    proteins=[]
    for protKey in pdbList:
        print protKey,reference,j
        if (protKey!=reference):
            curMol=B.PDBModel(protKey)
            curMol.renumberResidues()
            proteins.append(curMol)
            #curReducedMol=curMol.compress(curMol.maskProtein())
            rn= npy.unique(curMol["residue_number"])
            i=0
            for c,el in enumerate(seqDict[protKey]):
                if el!="-":
                    seqAlignment[j,c]=rn[i]
                    i+=1
                else :
                    seqAlignment[j,c]=-1
            j+=1

    return SequenceAlignment(reference,pdbList,seqAlignment,seqDict,submask,proteins)
#