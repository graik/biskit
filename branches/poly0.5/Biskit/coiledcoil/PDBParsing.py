from getcoil import createCandidatesFile

PDB_DATABASE_DIR = "/media/disk/PDB_Database/fs-sb/databases/pdb/pdbfiles"
OUTPUTFILE = "/media/disk/PDB_Database/fs-sb/databases/pdb/database.db"

createCandidatesFile(OUTPUTFILE,PDB_DATABASE_DIR)