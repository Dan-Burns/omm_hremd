import pdb2pqr
import pdbfixer
from openmm.app import PDBFile
import os


# Tools for quick system setup
# clean a structure file and protonate to specific pH with pdb2pqr
# using a Class might be the best approach

def Class SystemFromPDB():
    '''
    provide a pdb id or path to file 
    and output an OpenMM System object
    '''
    def __init__(self, pdb_file,
                 clean_pdb=f"{pdb_file}_cleaned.pdb",
                 protonated_pdb=f"{pdb_file}_protonated.pdb",
                 pH=7,
                 
                 )
    self.pdb_file = pdb_file

    self.pdb = pdbfixer.PDBFixer(self.pdb_file)

    self.pdb.removeHeterogens(keepWater=False)

     
    with open(f'{self.clean_pdb}', 'w') as f:
        PDBFile.writeFile(self.pdb.topology, self.pdb.positions, f)

    cmd = f"pdb2pqr30 --with-ph {self.pH} --ff={self.forcefield.upper()} {self.clean_pdb} {self.protonated_pdb}"

    #os.system(cmd)

    self.returned_output = subprocess.check_output(self.cmd, shell=True)