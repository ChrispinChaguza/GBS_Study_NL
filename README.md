## Data associated the the manuscript "Population genomics of Group B Streptococcus reveals the genetics of neonatal disease onset and meningeal invasion"

## Output files from the GWAS analysis using FaST-LMM
FaST-LMM output filename | Genetic variant type | Phenotype | Phenotype type
-- | -- | -- | -- 
GENE_final.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) 
GENE_final_GWAS.CNS_infection.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Meningeal (CNS) infection: Blood vs. CSF | Categorical (transformed) 
GENE_final_GWAS.ONSET_categorical.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Disease onset time (0-6 vs 7-89 days) | Categorical (transformed) |
GENE_final_GWAS.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) 
SNP_final.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | SNP | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) 
SNP_final_GWAS.CNS_infection.FaST-LMM.results.tsv.tar.gz | SNP | Meningeal (CNS) infection: Blood vs. CSF | Categorical 
SNP_final_GWAS.ONSET_categorical.FaST-LMM.results.tsv.tar.gz | SNP | Disease onset time (0-6 vs 7-89 days) | Categorical (transformed) 
SNP_final_GWAS.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | SNP | Disease onset time | Continuous (transformed) 
Unitigs_final.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) 
Unitigs_final_GWAS.CNS_infection.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Meningeal (CNS) infection: Blood vs. CSF | Categorical (transformed) 
Unitigs_final_GWAS.ONSET_categorical.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Disease onset time (0-6 vs 7-89 days) | Categorical (transformed) 
Unitigs_final_GWAS.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) 

## SNPs, accessory gene sequences, and unitig sequences
Variants/sequences | Filename | Description
-- | -- | -- 
Unitig presence and absence | Unitigs_presence_absence_matrix_all.tsv.tar.gz | All unitigs
Unitig presence and absence | Unitigs_presence_absence_matrix_filtered.tsv.tar.gz | Unitigs present in 5-95% isolates
Accessory gene presence and absence | accessory_gene_sequences.fa.tar.gz | All gene sequences
SNPs | SNP.vcf.tar.gz | All SNPs
Unitig sequences | Unitig_sequences.fa.tar.gz | All unitig sequences

## Scripts used to annotate SNPs and unitig sequences
Script name | Description
-- | -- 
annotate_SNPs.py | Generates a summary of gene features in a reference genome given SNP position
blast_annotate_fasta.py | Generates a summary of genetic features in GenBank-formatted reference genome(s) associated with given unitig sequences
GWAS_transform_quantitative_phenotype_RINT.R | R code for the rank-based inverse normal transformation of quantitative phenotypes/traits in GWAS (Original script written by Yuxuan Wang at Boston University; https://yuxuanstat.com/posts/2020/06/rank-based-inverse-normal-transformation/)
