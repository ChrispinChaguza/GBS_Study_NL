## Data associated the the manuscript "GWAS of Group B Streptococcus identifies genetic variation influencing neonatal disease onset and meningeal invasion"

## Output files from the GWAS analysis using FaST-LMM
FaST-LMM output filename | Genetic variant type | Phenotype | Phenotype type | Description
-- | -- | -- | -- | --
GENE_final.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) |
GENE_final_GWAS.CNS_infection.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Meningeal (CNS) infection: Blood vs. CSF | Categorical (transformed) |
GENE_final_GWAS.ONSET_categorical.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Disease onset time (0-6 vs 7-89 days) | Categorical (transformed) |
GENE_final_GWAS.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Gene presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) |
SNP_final.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | SNP | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) |
SNP_final_GWAS.CNS_infection.FaST-LMM.results.tsv.tar.gz | SNP | Meningeal (CNS) infection: Blood vs. CSF | Categorical |
SNP_final_GWAS.ONSET_categorical.FaST-LMM.results.tsv.tar.gz | SNP | Disease onset time (0-6 vs 7-89 days) | Categorical (transformed) |
SNP_final_GWAS.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | SNP | Disease onset time | Continuous (transformed) |
Unitigs_final.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) |
Unitigs_final_GWAS.CNS_infection.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Meningeal (CNS) infection: Blood vs. CSF | Categorical (transformed) |
Unitigs_final_GWAS.ONSET_categorical.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Disease onset time (0-6 vs 7-89 days) | Categorical (transformed) |
Unitigs_final_GWAS.ONSET_continuous.FaST-LMM.results.tsv.tar.gz | Unitigs presence/absence | Disease onset time (days from birth to GBS disease onset) | Continuous (transformed) |

## SNPs, accessory gene sequences, and unitig sequences
Variants/sequences | Filename | Description
-- | -- | -- 
Unitig presence and absence | Unitigs_presence_absence_matrix_all.tsv.tar.gz | All unitigs
Unitig presence and absence | Unitigs_presence_absence_matrix_filtered.tsv.tar.gz | Unitigs present in 5-95% isolates
Accessory gene presence and absence | accessory_gene_sequences.fa.tar.gz | All gene sequences
SNPs | SNP.vcf.tar.gz | All SNPs

## Scripts used to annotate SNPs and unitig sequences
Script name | Description
annotate_SNPs.py | Generates a summary of gene features in a reference genome given SNP position
blast_annotate_fasta.py | Generates a summary of genetic features in GenBank-formatted reference genome(s) associated with given unitig sequences
