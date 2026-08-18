# RIP-seq processing and analysis

This repository contains the processing and analysis workflow used for the
*Cryptosporidium parvum* F-RBP RIP-seq experiment. It starts from single-end
Ultima Genomics FASTQ files, performs trimming and STAR alignment, removes
reads likely to originate from unintended overlapping loci, quantifies the
remaining uniquely mapped reads with Salmon, and generates the downstream
analysis and coverage plots.

Except for the two RStudio sections identified below, run all commands from the
top-level repository directory.

## Requirements

Command-line tools:

- [FastQC](https://github.com/s-andrews/fastqc)
- [Cutadapt](https://github.com/marcelm/cutadapt)
- [gffread](https://github.com/gpertea/gffread)
- [STAR](https://github.com/alexdobin/STAR) 2.7.11b or later
- [samtools](https://github.com/samtools/samtools)
- [Salmon](https://github.com/COMBINE-lab/salmon)
- [deepTools](https://github.com/deeptools/deepTools)
- Python 3.9 or later; no non-standard Python packages are required

R packages used by `Analysis/RIP-seq_Analysis.R`:

- `ggplot2`
- `mixtools`
- `tximport`
- `svglite`

R/Bioconductor packages used by the coverage-plot scripts:

- `Gviz`
- `rtracklayer`
- `GenomicRanges`
- `IRanges`
- `svglite`

## Input files

The repository includes the CpBGF V1/V1.1 reference files obtained from the
[CpBGF repository](https://github.com/jkissing/CpBGF_Repository):

- `CpBGF_genome_v1.fasta`
- `CpBGF_genome_V1.1.gff`

Download the raw RIP-seq FASTQ files from NCBI BioProject `PRJNA1499133` and
place them in `RawData/`. The expected sample prefixes are:

- `1_WT_Bunchgrass_Input`
- `2_WT_Bunchgrass_IP`
- `3_F-RBP_Input`
- `4_F-RBP_IP`

Large sequencing files and generated outputs are excluded from Git by
`.gitignore`.

## 1. Trim and quality-check reads

These Ultima Genomics libraries were produced by PCR conversion of Illumina
libraries and contain excess sequence that must be removed. The command below
runs four samples concurrently, with six threads assigned to each Cutadapt and
FastQC process (up to 24 threads total).

```bash
mkdir -p trimmed_500

nohup bash -lc '
mkdir -p trimmed_500/qc_trimmed trimmed_500/cutadapt_logs

find ./RawData -maxdepth 1 -name "*.fastq.gz" -print0 \
| xargs -0 -P 4 -I {} bash -lc "
  f=\"{}\"
  base=\$(basename \"\$f\" .fastq.gz)

  cutadapt -j 6 \\
    -u 22 \\
    -a NNNNNNNNNNNNNNAGATCGGAAGAG \\
    -q 20 \\
    -m 50 \\
    -l 500 \\
    -o \"trimmed_500/\${base}.trim.fastq.gz\" \"\$f\" \\
    > \"trimmed_500/cutadapt_logs/\${base}.log\" 2>&1

  fastqc -t 6 \\
    -o trimmed_500/qc_trimmed \\
    \"trimmed_500/\${base}.trim.fastq.gz\"
"
' > trimmed_500/trim_500_nohup.log 2>&1 &
```

## 2. Prepare the alignment reference

Create a GFF containing only protein-coding genes and lncRNAs:

```bash
GFF="./CpBGF_genome_V1.1.gff"
OUT="./CpBGF_genome_V1.1_mRNA_lncRNA_only.gff"

awk -F '\t' '
BEGIN { OFS="\t" }
$0 ~ /^#/ { next }
{
  split($9, fields, ";")
  delete attr
  for (i in fields) {
    split(fields[i], key_value, "=")
    attr[key_value[1]] = key_value[2]
  }

  if ($3 == "gene" &&
      (attr["gene_biotype"] == "protein_coding" ||
       attr["gene_biotype"] == "lncRNA")) {
    keep_gene[attr["ID"]] = 1
    print
    next
  }

  parent = attr["Parent"]
  if (parent in keep_gene) {
    keep_transcript[attr["ID"]] = 1
    print
    next
  }

  if (parent in keep_transcript) print
}
' "$GFF" > "$OUT"
```

Remove seven lncRNAs that directly overlap antisense rRNA loci. This version
uses only POSIX-compatible awk features.

```bash
awk -F '\t' '
BEGIN {
  OFS="\t"
  split("cpbgf_1003875 cpbgf_7005560 cpbgf_7005570 cpbgf_7005840 cpbgf_7005870 cpbgf_8005465 cpbgf_8005505", loci, " ")
  for (i in loci) drop[loci[i]] = 1
}
{
  locus = ""
  n = split($9, fields, ";")
  for (i = 1; i <= n; i++) {
    split(fields[i], key_value, "=")
    if (key_value[1] == "locus_tag") locus = key_value[2]
  }
  if (!(locus in drop)) print
}
' ./CpBGF_genome_V1.1_mRNA_lncRNA_only.gff \
  > ./CpBGF_genome_V1.1_mRNA_lncRNA_only_clean.gff
```

Convert the cleaned GFF to GTF and generate the STAR index:

```bash
gffread ./CpBGF_genome_V1.1_mRNA_lncRNA_only_clean.gff \
  -T \
  -o ./CpBGF_genome_V1.1_mRNA_lncRNA_only_clean.gtf

mkdir -p ./STAR_mRNA_lncRNA_only_clean_500bp

STAR --runThreadN 32 \
  --runMode genomeGenerate \
  --genomeDir ./STAR_mRNA_lncRNA_only_clean_500bp/ \
  --genomeFastaFiles ./CpBGF_genome_v1.fasta \
  --sjdbGTFfile ./CpBGF_genome_V1.1_mRNA_lncRNA_only_clean.gtf \
  --sjdbOverhang 499 \
  --genomeSAindexNbases 10
```

## 3. Align reads with STAR

`STAR_alignment_script.sh` defaults to four concurrent STAR jobs with eight
threads per job. Its `--limitBAMsortRAM 100000000000` setting permits up to
100 GB of BAM-sorting RAM per STAR process, or up to 400 GB across four
concurrent jobs. Adjust `JOBS`, `THREADS_PER_JOB`, `STAR_INDEX`, `OUTDIR`, and
the RAM limit before running on a smaller system.

The script uses `--quantTranscriptomeSAMoutput BanSingleEnd`, which allows
indels and soft clips in the transcriptome BAM. Change this to
`BanSingleEnd_ExtendSoftclip` only if soft-clipped sequence should be extended
in the transcriptome projection.

```bash
chmod +x ./STAR_alignment_script.sh
nohup bash -lc './STAR_alignment_script.sh' \
  > STAR_Alignment_Log.log 2>&1 &
```

STAR writes its outputs to `star_out_mRNA_lncRNA_only_clean_500bp/`, including
paired coordinate and transcriptome BAMs for each sample.

## 4. Remove likely misaligned reads

The cleanup wrapper uses only standard Python plus samtools. Its default paths
already match this repository layout:

- full annotation: `CpBGF_genome_V1.1.gff`
- alignment annotation: `CpBGF_genome_V1.1_mRNA_lncRNA_only_clean.gff`
- STAR BAM directory: `star_out_mRNA_lncRNA_only_clean_500bp/`
- output directory: `misaligned_read_filter_500bp/`

The cleanup removes all transcriptome alignments belonging to reads selected by
the rules outlined below. It also writes coordinate BAMs containing the
discarded reads for inspection.

Rules identifying alignments to be removed:

1. **Always-remove loci:** remove reads aligned to any locus listed in
   `ALWAYS_REMOVE_LOCI` (by default, `cpbgf_3001410`). This rule is unconditional
   and does not depend on read orientation.
2. **Special-locus overlap:** remove reads whose projected genomic alignment
   overlaps `SPECIAL_LOCUS` (by default, `cpbgf_300675`). This direct-locus rule
   is also unconditional with respect to read orientation.
3. **Same-orientation unintended overlap:** remove reads when all of the
   following are true:
   - the read is aligned to a transcript in the restricted reference GFF;
   - that reference gene overlaps an unintended gene on the same strand in the
     full GFF;
   - the read's projected genomic alignment lies within the unintended gene,
     allowing `EDGE_TOLERANCE_NT` bases at each edge (5 nt by default); and
   - the read orientation matches `SAME_ORIENTATION_READ_STRAND` (`reverse` by
     default for Salmon library type `SR`; use `forward` for `SF`).

The final removal set is the union of these rules, so a read triggering more
than one rule is removed only once. From a biological standpoint, these selection criteria are designed to identify reads that STAR mapped to an mRNA or lncRNA, but that most likely actually originated from a rRNA locus that is nested within said mRNA or lncRNA, in the same orientation. By visual inspection of alignments using IGV, we found several of those cases that lead to substantial distortion of the final TPM estimates across samples.

```bash
chmod +x ./filter_misaligned_reads.sh ./filter_misaligned_reads.py
bash ./filter_misaligned_reads.sh
```

The transcriptome BAMs used by the next step are:

```text
misaligned_read_filter_500bp/*.Aligned.toTranscriptome.clean.out.bam
```

## 5. Retain only genomically unique reads

This step identifies reads whose primary genomic STAR alignment has `NH:i:1`,
then retains all transcript projections for those reads so Salmon can resolve
compatible transcript assignments.

```bash
chmod +x ./make_clean_unique_transcriptome_bams.sh
bash ./make_clean_unique_transcriptome_bams.sh
```

Outputs are written to:

```text
unique_transcriptome_bams_500bp/*.Aligned.toTranscriptome.clean.unique.out.bam
```

## 6. Quantify transcripts with Salmon

Generate the full V1.1 transcriptome FASTA. Using the full GFF here is
intentional: transcripts excluded from the STAR reference remain represented
as zero-alignment targets rather than disappearing from the quantification
reference.

```bash
gffread ./CpBGF_genome_V1.1.gff \
  -g ./CpBGF_genome_v1.fasta \
  -w ./CpBGF_genome_V1.1.transcripts.fasta
```

Run Salmon in alignment-based mode. `-l SR` describes a strand-specific
single-end library whose reads align in reverse orientation relative to their
transcripts; `-s` enables sequence-bias correction.

```bash
mkdir -p ./salmon_quant ./salmon_logs

find ./unique_transcriptome_bams_500bp \
  -maxdepth 1 \
  -name '*.Aligned.toTranscriptome.clean.unique.out.bam' \
  -print0 \
| xargs -0 -P 4 -I {} bash -lc '
    bam="$1"
    filename=$(basename "$bam")
    base=${filename%.Aligned.toTranscriptome.clean.unique.out.bam}

    salmon quant \
      -t ./CpBGF_genome_V1.1.transcripts.fasta \
      -l SR \
      -s \
      -a "$bam" \
      -o "./salmon_quant/${base}" \
      > "./salmon_logs/salmon_${base}.log" 2>&1
  ' _ {}
```

## 7. Transcript-to-gene mapping

`Analysis/tx2gene.csv` is included in the repository and covers every
transcript in the full V1.1 transcriptome FASTA. It can be regenerated with:

```bash
awk -F '\t' '
BEGIN {
  OFS=","
  print "transcript_id", "gene_id"
}
$3 == "mRNA" || $3 == "lnc_RNA" || $3 == "rRNA" ||
$3 == "tRNA" || $3 == "snRNA" || $3 == "snoRNA" ||
$3 == "SRP_RNA" {
  transcript = ""
  gene = ""
  n = split($9, fields, ";")
  for (i = 1; i <= n; i++) {
    split(fields[i], key_value, "=")
    if (key_value[1] == "ID") transcript = key_value[2]
    if (key_value[1] == "locus_tag") gene = key_value[2]
    if (key_value[1] == "Parent" && gene == "") {
      gene = key_value[2]
      sub(/^gene-/, "", gene)
    }
  }
  if (transcript != "" && gene != "") print transcript, gene
}
' ./CpBGF_genome_V1.1.gff > ./Analysis/tx2gene.csv
```

## 8. Run the R analysis

Open `Analysis/Analysis.Rproj` in RStudio, then open and run
`Analysis/RIP-seq_Analysis.R`. Opening the project is important: it makes
`Analysis/` the working directory, so the script can find `../salmon_quant/`
and the CSV files stored beside the R script.

The script checks for all four expected Salmon outputs before analysis. It
generates Figure 5B, Supplementary Figure 5, and the data underlying
Supplementary Table 2 in `Analysis/output/`.

## 9. Generate coordinate BAMs and BigWig coverage

Create coordinate BAMs containing exactly the genomically unique reads retained
in the final transcriptome BAMs:

```bash
chmod +x ./make_clean_unique_coordinate_bams.sh
bash ./make_clean_unique_coordinate_bams.sh
```

Generate CPM-normalized BigWig files. The settings below run four samples in
parallel with 16 threads per sample.

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

INDIR="./coord_bams_500_clean_unique"
JOBS=4
THREADS_PER_JOB=16
export THREADS_PER_JOB

bams=("$INDIR"/*Aligned.sortedByCoord.clean.unique.out.bam)
if (( ${#bams[@]} == 0 )); then
  echo "ERROR: No matching BAM files found in $INDIR" >&2
  exit 1
fi

printf '%s\0' "${bams[@]}" \
| xargs -0 -P "$JOBS" -I {} bash -lc '
    set -euo pipefail

    bam="$1"
    output_dir="$2"
    filename=$(basename "$bam")
    base=${filename%.Aligned.sortedByCoord.clean.unique.out.bam}
    base=${base%.trim}
    out="${output_dir}/${base}.bw"

    echo "[bamCoverage] $base"
    bamCoverage \
      -b "$bam" \
      -o "$out" \
      -p "$THREADS_PER_JOB" \
      --normalizeUsing CPM \
      --binSize 1 \
      --smoothLength 50
  ' _ {} "$INDIR"
```

Generate the full GTF used for gene annotation in the coverage plots:

```bash
gffread ./CpBGF_genome_V1.1.gff \
  -T \
  -o ./CpBGF_genome_V1.1.gtf
```

Finally, open `Coverage_Plots/CoveragePlots.Rproj` in RStudio and run:

- `CoveragePlot_cgd4_1910.R`
- `CoveragePlot_cgd6_2450.R`

The project sets `Coverage_Plots/` as the working directory, allowing the
scripts to find the GTF and BigWig files one directory above it.
