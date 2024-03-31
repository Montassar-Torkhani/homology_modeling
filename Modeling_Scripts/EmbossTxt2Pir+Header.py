from Bio import AlignIO

def get_structure_type(pdb_file):
    """Extracts the structure type from the PDB file."""
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("EXPDTA"):
                if "X-RAY DIFFRACTION" in line:
                    return "structureX"
                elif "SOLUTION NMR" in line:
                    return "structureN"
    return None

def get_structure_id(pdb_file):
    """Extracts the structure ID from the PDB file."""
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("HEADER"):
                return line.rstrip().split()[-1]
    return None

def get_protein_info(pdb_file):
    """Extracts protein name and organism from the PDB file."""
    protein_name = None
    organism = None
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("COMPND"):
                parts = line.split()
                if len(parts) >= 4 and parts[1].isdigit() and parts[2] == "MOLECULE:":
                    name = parts[-1].strip(";")
                    if len(name) >= 2 and name[0] in ('"', "'", "(") and name[-1] in ('"', "'", ")"):
                        protein_name = name[1:-1]  # Remove first and last characters
                    else:
                        protein_name = name
                                        # Check if the last word is 'A' or 'B', if so, take the second last word as the protein name
                    if protein_name and protein_name[-1] in ('A', 'B'):
                        
                        protein_name = parts[-2]
            elif line.startswith("SOURCE"):
                parts = line.split()
                if len(parts) >= 5 and parts[1].isdigit() and parts[2] == "ORGANISM_SCIENTIFIC:":
                    organism = " ".join(parts[3:]).strip(";")
            if protein_name and organism:
                break
    return protein_name, organism

def get_resolution_and_rfactor(pdb_file):
    """Extracts resolution and R-factor from the PDB file."""
    resolution = None
    rfactor = None
    with open(pdb_file, "r") as f:
        for line in f:
            if "REMARK" in line:
                if "RESOLUTION RANGE HIGH" in line:
                    resolution = line.split(":")[1].strip()
                elif "R VALUE" in line and "WORKING SET" in line:
                    rfactor = line.split(":")[1].strip()
            if resolution and rfactor:
                break
    return resolution, rfactor

def get_first_and_last_atom_info(pdb_file):
    first_residue_position = None
    first_residue_chain = None
    last_residue_position = None
    last_residue_chain = None
    
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                if not first_residue_position:
                    first_residue_position = int(parts[5])
                    first_residue_chain = parts[4]
            elif line.startswith("TER"):
                parts = line.split()
                last_residue_chain = parts[3]
                last_residue_position = int(parts[4])
                

    return first_residue_position, first_residue_chain, last_residue_position, last_residue_chain

def emboss_to_pir(emboss_file, pdb_file, pir_file):
    """Converts an EMBOSS alignment file to a custom PIR-like format."""
    
    alignment = AlignIO.read(emboss_file, "emboss")
    
    structure_type = get_structure_type(pdb_file)
    structure_id = get_structure_id(pdb_file)
    protein_name, organism = get_protein_info(pdb_file)
    resolution, rfactor = get_resolution_and_rfactor(pdb_file)
    first_residue_position, first_residue_chain, last_residue_position, last_residue_chain = get_first_and_last_atom_info(pdb_file)
    
    with open(pir_file, "w") as f:
        if structure_type == "structureX":
            f.write(f">P1;{structure_id}_clean\n")
            f.write(f"{structure_type}:{structure_id}_clean:{first_residue_position}    :{first_residue_chain}:{last_residue_position}  :{last_residue_chain}:{protein_name}:{organism}: {resolution}: {rfactor}\n")
        elif structure_type == "structureN":
            f.write(f">P1;{structure_id}\n")
            f.write(f"{structure_type}:{structure_id}:{first_residue_position}    :{first_residue_chain}:{last_residue_position}  :{last_residue_chain}:{protein_name}:{organism}::\n")

        for i, record in enumerate(alignment, start=1):
            if i == 1:
                sequence = str(record.seq).replace("-", "-")  # Get sequence for the first record
                f.write(f"{sequence}*\n")
                f.write("\n") 
            else:
                f.write(f">P1;{record.id}\n")
                f.write(f"sequence:{record.id}::::::::\n")  # Custom header format for sequence 2 
                sequence = str(record.seq).replace("-", "-")  # Replace gaps with dashes
                f.write(f"{sequence}*\n")
                
# Replace with your actual input and output file paths
pdb_file = '1c8q.pdb'
input_txt_alignment = 'emboss_pairwise_alignment5.txt'
pir_file = 'alignment.pir'

emboss_to_pir(input_txt_alignment, pdb_file, pir_file)
