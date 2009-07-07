from getcoil import createCandidatesFile

PDB_DATABASE_DIR = "/fs-sb/databases/pdb/pdbfiles"
OUTPUTFILE = "/home/nisusers/vgil/coils.db"

createCandidatesFile(OUTPUTFILE,PDB_DATABASE_DIR)