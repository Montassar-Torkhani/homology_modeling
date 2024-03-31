import os,re

# Define paths for PDB and FASTA files
pdb_file = '3blk.pdb'
fasta_file = '3blk_Fa.fasta'

aa3to1 = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
    'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
    'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
    'MSE': 'M',
}

ca_pattern = re.compile(r"^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
filename = os.path.basename(pdb_file).split('.')[0]
chain_dict = dict()
chain_list = []

with open(pdb_file, 'r') as fp:
    for line in fp.read().splitlines():
        if line.startswith("ENDMDL"):
            break
        match_list = ca_pattern.findall(line)
        if match_list:
            resn = match_list[0][0] + match_list[0][2]
            chain = match_list[0][1] + match_list[0][3]
            if chain in chain_dict:
                chain_dict[chain] += aa3to1[resn]
            else:
                chain_dict[chain] = aa3to1[resn]
                chain_list.append(chain)

with open(fasta_file, 'w') as fasta_fp:
    for chain in chain_list:
        fasta_fp.write('>%s:%s\n%s\n' % (filename, chain, chain_dict[chain]))

print(f"FASTA sequence extracted and saved to {fasta_file}")
