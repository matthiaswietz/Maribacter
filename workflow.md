## CAZyme & growth data

Analyses and figures from the paper can be reproduced using this [Rscript](./code/Maribacter.R) and [input files](./data/Rstats) 

## Pangenome analysis

- We use protein-fasta files (faa--translated cds) from [23 strains](./data/pangenome_faa)

- For processing, transfer these to a folder "Maribacter" inside your Linux working directory.  

- Then run OrthoFinder v2.3.3 as follows:  

`OrthoFinder-2.3.3/orthofinder -o ~/Maribacter/results -f ~/Maribacter/ -M msa -t 17`

- The output file `Orthogroups.txt` is then modified as follows:

`cd ~/Maribacter/results/Results/Orthogroups`  
`grep "OG*" Orthogroups.txt | sed 's/://' > Orthogroups_2.txt`

####################################################

- To derive core / pan genes, run this [perl script](./code/Pangenome.sh) from your working directory    
(original script available from https://zenodo.org/record/1010076/files/PanGenome_analysis_host_specific.tar.gz)

- We define three groups in [this file](./data/Strain_list.txt):
   - Maribacter_621: Group1 (6 strains)
   - Maribacter_count: Group2 (9 strains)
   - Zobellia_count: Group3 (8 strains)

- Then execute:   
`Pangenome.sh ~/Maribacter/Strain_list.txt ~/Maribacter/results/Results/Orthogroups/Orthogroups_2.txt`

- This creates a file with core/pan gene lists. We furthermore extract the strains-specific unique genes (assigned to Orthogroups when detected in >1 copy in a given strain).

- Unassigned unique genes are written into
`~/Maribacter/results/Results/Orthogroups/Orthogroups_UnassignedGenes.tsv`

- This is formatted as follows:
`grep "OG*" Orthogroups_UnassignedGenes.tsv | sed 's/\t\+/\t/g' | sed '1d' > Orthogroups_UnassignedGenes_LocusTag.txt`

- Unique genes are derived from `Orthogroups.GeneCount.tsv` by executing  
`head -1 Orthogroups.GeneCount.tsv | sed 's/[^\t]//g' | wc -c`
`for i in {2.."n-1"}
do
awk -v i="$i" '$i== $"n" { print $0}' Orthogroups.GeneCount.tsv >> UniqueOrthogroups.tsv
done

The output is then combined in Rstudio with dbCAN2 or KEGG/KAAS assignments by merging with the locus_tag.
