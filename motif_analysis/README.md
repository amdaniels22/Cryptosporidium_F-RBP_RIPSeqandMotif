# Y-BOX Binding Motif Analysis Pipeline

A comprehensive bioinformatics pipeline for identifying and analyzing Y-BOX binding motifs in *C. parvum* (CpBGF) mRNA sequences using MEME suite tools and R-based visualization.

## Overview

This pipeline performs motif discovery and enrichment analysis on RIP (RNA Immunoprecipitation) experimental data. It extracts coding sequences from annotated genomes, runs MEME motif discovery, performs FIMO scanning against the whole transcriptome, and visualizes results using R.

## Prerequisites

### System Requirements
- Linux/Unix environment
- R version ≥ 4.5.1
- Python 3.x with pandas and numpy

### Software Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| AGAT | 1.7.0 | `conda install -c bioconda agat` |
| MEME Suite | 5.5.9 | `conda install -c bioconda meme` |
| Samtools | Latest | `conda install -c bioconda samtools` |
| R packages | - | See R section below |

### R Packages
```r
# Install if missing
install.packages(c("Seurat", "ggplot2", "ggrepel", "dplyr", 
                   "pheatmap", "scales", "FactoMineR", "factoextra"))
```

## Directory Structure

```
.
├── data/
│   ├── CpBGF_genome_V1.gff          # Reference genome annotation
│   ├── CpBGFT2T.fasta                # Reference genome sequences
│   └── rip_list.txt                  # Input gene names from RIP experiment
├── scripts/
│   ├── make_master_matrix.py         # Python script for matrix generation
│   └── motif_analysis_script.R       # R script for visualization
├── output/
│   ├── cds_rna_motif_1_fimo/         # FIMO output directory
│   ├── rip_list_mrna_cds_rna_motif/  # MEME output directory
│   └── ...                            # Other intermediate files
└── README.md                          # This file
```

## Workflow

### Step 1: Filter Genes from RIP List

Extract genes of interest from the genome annotation:

```bash
agat_sp_filter_feature_from_keep_list.pl \
    --gff CpBGF_genome_V1.gff \
    --keep_list rip_list.txt \
    -o rip_list.gff
```

### Step 2: Generate Statistics

```bash
# General GFF statistics
agat_sp_statistics.pl \
    --gff rip_list.gff \
    --gs 9259183 \
    -d --output rip_list_gff_report.txt

# Functional statistics
agat_sp_functional_statistics.pl \
    --gff rip_list.gff \
    --gs 9259183
```

### Step 3: Separate mRNA from lncRNA

```bash
agat_sp_separate_by_record_type.pl \
    --gff rip_list.gff
```

*Note: This creates `split_result/mrna.gff` for downstream analysis.*

### Step 4: Extract Sequences

Extract various genomic regions for motif analysis:

```bash
# CDS sequences
agat_sp_extract_sequences.pl \
    -g split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t cds -o rip_list_mrna_cds.fasta

# CDS protein sequences
agat_sp_extract_sequences.pl \
    -g split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t cds -p -o rip_list_mrna_prot.fasta

# 5' UTR
agat_sp_extract_sequences.pl \
    -g split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t five_prime_utr --merge -o rip_list_mrna_5utr.fasta

# 3' UTR
agat_sp_extract_sequences.pl \
    -g split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t three_prime_utr --merge -o rip_list_mrna_3utr.fasta

# 150bp upstream region
agat_sp_extract_sequences.pl \
    -g split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t gene -5 150 --full -o rip_list_mrna_5up.fasta

# 150bp downstream region
agat_sp_extract_sequences.pl \
    -g split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t gene -3 150 --full -o rip_list_mrna_3down.fasta
```

### Step 5: Run MEME Motif Discovery

Run MEME on each extracted sequence type:

```bash
# 3' downstream (DNA)
meme rip_list_mrna_3down.fasta \
     -mod anr -dna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -revcomp -oc rip_list_mrna_3down_motif

# 5' UTR (DNA)
meme rip_list_mrna_5utr.fasta \
     -mod anr -dna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -revcomp -oc rip_list_mrna_5utr_dna_motif

# 3' UTR (DNA)
meme rip_list_mrna_3utr.fasta \
     -mod anr -dna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -revcomp -oc rip_list_mrna_3utr_dna_motif

# CDS (DNA)
meme rip_list_mrna_cds.fasta \
     -mod anr -dna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -revcomp -oc rip_list_mrna_cds_dna_motif

# Protein sequences
meme rip_list_mrna_prot.fasta \
     -mod anr -protein -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -oc rip_list_mrna_prot_motif

# CDS (RNA)
meme rip_list_mrna_cds.fasta \
     -mod anr -rna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -oc rip_list_mrna_cds_rna_motif

# 3' UTR (RNA)
meme rip_list_mrna_3utr.fasta \
     -mod anr -rna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -oc rip_list_mrna_3utr_rna_motif

# 5' UTR (RNA)
meme rip_list_mrna_5utr.fasta \
     -mod anr -rna -nmotifs 5 -minw 4 -maxw 9 -maxsize 999999999 \
     -oc rip_list_mrna_5utr_rna_motif
```

**Note:** After inspection, only `rip_list_mrna_cds_rna_motif` typically produces motifs passing the default p-value cutoff.

### Step 6: Motif Comparison with Tomtom

Compare discovered motifs against known databases:

```bash
tomtom \
    -no-ssc -oc rip_list_mrna_cds_rna_motif/tomtom \
    -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 \
    cds_rna_motif_1.meme \
    db/JASPAR/JASPAR2026_CORE_vertebrates_non-redundant.meme \
    db/EUKARYOTE/jolma2013.meme \
    db/MOUSE/uniprobe_mouse.meme
```

### Step 7: Genome-wide FIMO Scanning

Prepare for FIMO scanning of the entire transcriptome:

```bash
# Download reference from: https://github.com/jkissing/CpBGF_Repository

# Split genome GFF
agat_sp_separate_by_record_type.pl \
    -o fullgff_split_result --gff CpBGF_genome_V1.gff

# Extract all CDS sequences
agat_sp_extract_sequences.pl \
    -g fullgff_split_result/mrna.gff \
    -f CpBGFT2T.fasta \
    -t cds -o all_mrna_cds.fasta

# Index FASTA file
samtools faidx all_mrna_cds.fasta

# Create length file for Python script
cut -f 1,2 all_mrna_cds.fasta.fai \
    | sed 's/rna-gnl|cpbgf|//g' > all_mrna_cds_length.txt
```

Run FIMO:

```bash
fimo -oc cds_rna_motif_1_fimo \
    --verbosity 5 \
    cds_rna_motif_1.meme \
    all_mrna_cds.fasta
```

### Step 8: Process FIMO Output

Generate gene count tables:

```bash
# Count motif occurrences per gene
cat fimo.tsv | tail -n+2 \
    | cut -f 3 | sort -V \
    | sed 's/rna-gnl|cpbgf|//g' | uniq -c \
    | sort -k1,1rn | awk '$2!~"^#"' | awk '$2!=""' \
    | awk '{print $2}' | awk -F '-' '{print $1}' \
    | sort -V | uniq -c \
    | sort -k1,1rn | awk '{print $2}' | sort -V
```

### Step 9: Generate Master Matrix

Prepare input files for Python processing:

```bash
# Gene lengths
cut -f 1,2 all_mrna_cds.fasta.fai \
    | sed 's/rna-gnl|cpbgf|//g' > all_mrna_cds_length.txt

# Enriched gene counts
cat fimo.tsv | tail -n+2 \
    | cut -f 3 | sort -V \
    | sed 's/rna-gnl|cpbgf|//g' | uniq -c \
    | sort -k1,1rn | awk '$2!~"^#"' | awk '$2!=""' \
    | awk '{print $2,$1}' OFS='\t' > enriched_genes_counts.txt

# Input gene names
cat rip_list.txt \
    | sed 's/gene-//g' | awk '{print $1"-RA"}' > input_gene_name.txt

# Clean FIMO TSV
awk -F'\t' -vOFS='\t' '!/^#/ && NF {sub(/.*\|/,"",$3); print}' fimo.tsv > fimo_cleaned.tsv
```

Run Python script to create master matrix:

```bash
python3 make_master_matrix.py
```

### Step 10: Generate Subsets for R Analysis

Create gene lists from master matrix:

```bash
# All genes with motifs
awk 'NR>1 {print $1}' master_data_matrix.tsv > fimo_iowa2_name.txt

# Gene-motif count pairs for pie charts
awk 'NR>1 {print $1"\t"$2}' master_data_matrix.tsv > piedata.tsv

# Subset by motif count
awk 'NR>1 && $2==1 {print $1}' master_data_matrix.tsv > fimo_1motif_iowa2_name.txt
awk 'NR>1 && $2==2 {print $1}' master_data_matrix.tsv > fimo_2motif_iowa2_name.txt
awk 'NR>1 && $2==3 {print $1}' master_data_matrix.tsv > fimo_3motif_iowa2_name.txt
awk 'NR>1 && $2>=4 {print $1}' master_data_matrix.tsv > fimo_4plusmotif_iowa2_name.txt
```

### Step 11: R Visualization

Run the R analysis script in either terminal R, or interactively via RStudio:

```bash
Rscript motif_analysis_script.R
```

This generates:
- UMAP feature plots showing motif enrichment
- Module score density plots
- Pie charts of motif occurrence distribution

## Output Files

| File | Description |
|------|-------------|
| `master_data_matrix.tsv` | Main matrix with gene names and motif counts |
| `fimo_iowa2_name.txt` | All genes with detected motifs |
| `fimo_1motif_iowa2_name.txt` | Genes with exactly 1 motif |
| `fimo_2motif_iowa2_name.txt` | Genes with exactly 2 motifs |
| `fimo_3motif_iowa2_name.txt` | Genes with exactly 3 motifs |
| `fimo_4plusmotif_iowa2_name.txt` | Genes with 4+ motifs |
| `piedata.tsv` | Gene-motif count pairs for visualization |

## 🔗 References

- **Genome Reference**: [CpBGF_Repository](https://github.com/jkissing/CpBGF_Repository)
- **MEME Suite**: [meme-suite.org](https://meme-suite.org/)
- **JASPAR Database**: [jaspar.genereg.net](https://jaspar.genereg.net/)

## License

This project is for research purposes only.


