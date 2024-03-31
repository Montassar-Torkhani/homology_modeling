## Homology Modeling of Bacterial Extremophile Alpha-Amylase (1AQM)

This repository details the investigation of a bacterial extremophile alpha-amylase from Pseudoalteromonas haloplanktis (1AQM) using homology modeling with Modller. The project explores the effectiveness of identifying a suitable template protein despite a potentially low sequence identity between the target and templates.

### Methodology

Template Selection:

Four templates from diverse organisms were chosen for initial exploration:
1C8Q (Homo sapiens)
1VAH (Drosophila melanogaster)
3VM5 (Sus scrofa)
8ORP (Oryzias latipes)
By considering factors like sequence similarity and functional relevance, the most suitable template(s) were selected for further analysis.
Sequence Preparation:

PDB files of the chosen template(s) were converted to FASTA format using a custom Python script (clean_PDB.py). This ensures the sequences are in the appropriate format for downstream analysis.
Pairwise Alignment:

EMBOSS Needle was employed to perform pairwise alignments between the target sequence (1AQM) and the chosen template(s).
Another Python script (EmbossTxt2Pir+Header.py) was used to convert the alignment results from EMBOSS Needle into PIR format, potentially adjusting sequence headers to ensure proper formatting for homology modeling.
Homology Modeling:

Modller software was used to construct homology models of the target protein (1AQM) based on the selected template(s) and sequence alignment(s).
A dedicated Python script (Modeller_script.py) facilitated this process, potentially allowing for customization of modeling parameters.
Model Evaluation:

Various parameters were used to assess the quality of the generated models:
Ramachandran plot analysis identifies the distribution of amino acid backbone torsion angles, indicating regions of potentially favorable or unfavorable conformations.
QMEANDisCo Global score provides an overall estimate of the model's quality.
DOPE score and Z-score RMSD are additional metrics used to evaluate model quality.
Based on a combination of these metrics, the best models were selected for further investigation.

### Results and Discussion

Following model evaluation using Ramachandran plots, QMEAN scores, DOPE scores, and RMSD values, two models emerged as the most favorable:

1C8Q_Model: This model was built based on the template 1C8Q .
8ORP_Model: This model was built based on the template 8ORP .
Structural Comparison and Conservation:

To assess the accuracy of the models, we performed a comparison alignment. The results were:

1C8Q_Model vs. Target (1AQM): RMSD = 0.772
Template 1C8Q vs. Target (1AQM): RMSD = 0.659
While the template 1C8Q exhibited a slightly lower RMSD compared to the target protein (1AQM), both the 1C8Q_Model and 8ORP_Model achieved relatively low RMSD values, indicating good structural agreement with the target.

This success highlights the principle of conservation between protein sequence and structure. Even with potentially low sequence identity between the target and templates, key amino acids responsible for the overall protein fold can be conserved.  This allows for structural prediction using homology modeling despite sequence variations.


## Software Used

- Anaconda: Python distribution and package manager (https://www.anaconda.com/)
- Modller: Software for homology modeling (https://modbase.compbio.ucsf.edu/)
- Biopython : Python library for computational molecular biology (https://biopython.org/)
