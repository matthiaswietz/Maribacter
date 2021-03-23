## PANGENOME ANALYSIS 

We use protein-fasta files (faa--translated cds) from [23 strains](./data/pangenome_faa)

For processing, transfer these to a folder "Maribacter" inside your Linux working directory.  

Then run OrthoFinder v2.3.3 as follows:  

OrthoFinder-2.3.3/orthofinder -o ~/Maribacter/results -f ~/Maribacter/ -M msa -t 17

The resulting file ./Results/Orthogroups.txt, is then modified as follows:

grep "OG*" Orthogroups.txt | sed 's/://' > Orthogroups_2.txt

####################################################

To derive core / pan genes from the results, run this [perl script](./code/Pangenome.sh)     
(modified from https://zenodo.org/record/1010076#.YFKs_q_7Q2z)

For this, we define three groups in a file Strain_list.txt:
- Maribacter_621: Group1 (6 strains: 62-1 )
- Maribacter_count: Group2 and Group3 (9 strains)
- Zobellia_count: Outgroup (8 strains)

Then execute: 
Pangenome.sh ./Maribacter/Strain_list.txt ./Maribacter/Results/Orthogroups/Orthogroups_2.txt
