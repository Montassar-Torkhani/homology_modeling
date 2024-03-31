from Bio.PDB import PDBParser, Select, PDBIO

class CleanPDB(Select):
    def accept_residue(self, residue):
        # Remove HETATM records
        if residue.id[0].startswith('H_') or residue.id[0] == 'W':
            return 0
        # Remove hydrogen atoms
        elif any(atom.element == 'H' for atom in residue.get_atoms()):
            return 0
        else:
            return 1

    def accept_atom(self, atom):
        # Remove hydrogen atoms
        if atom.element == 'H':
            return 0
        else:
            return 1
        
pdb_parser = PDBParser()
structure = pdb_parser.get_structure('PDB_structure', '1c8q.pdb')

# Initialize PDBIO object
pdb_io = PDBIO()
pdb_io.set_structure(structure)
pdb_io.save('1c8q_clean.pdb', CleanPDB())

# Print a success message
print("PDB file  was successfully cleaned and saved ")
