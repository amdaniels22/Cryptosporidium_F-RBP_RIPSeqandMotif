#!/usr/bin/env bash
# Align each trimmed, single-end FASTQ file to the STAR genome index.
#
# This script is intentionally a small batch driver: it discovers FASTQs in
# ./trimmed_500, starts several independent STAR jobs, and optionally indexes
# each coordinate-sorted BAM when samtools is available.
set -euo pipefail

# Number of samples to process at once and STAR threads per sample.  The
# product (JOBS * THREADS_PER_JOB) is the approximate CPU demand.
JOBS=4
THREADS_PER_JOB=8
# This directory must already contain a STAR genome generated from the GFF and
# FASTA described in the README.
STAR_INDEX="./STAR_mRNA_lncRNA_only_clean_500bp/"
# STAR writes all per-sample output files beneath this directory.
OUTDIR="star_out_mRNA_lncRNA_only_clean_500bp"

mkdir -p "${OUTDIR}"

# mapfile collects the sorted list into a Bash array, preserving filenames as
# individual arguments even if they contain unusual characters.
mapfile -t FASTQS < <(
  find ./trimmed_500 -maxdepth 1 -name "*.fastq.gz" -print | sort
)

if [[ ${#FASTQS[@]} -eq 0 ]]; then
  echo "No SE fastq files found in ./trimmed_500."
  exit 1
fi

# xargs supplies one FASTQ to each worker.  The nested shell is used so that
# each worker has its own strict-error settings and local variables.
printf '%s\n' "${FASTQS[@]}" \
| xargs -P "${JOBS}" -I{} bash -lc '
  set -euo pipefail

  # The final "_ {}" below makes {} arrive as $1 in this nested shell.
  FQ="$1"
  BASE=$(basename "$FQ" .fastq.gz)

  # STAR appends names such as Aligned.sortedByCoord.out.bam to this prefix.
  PREFIX="'"${OUTDIR}"'/${BASE}."
  echo "[STAR] ${BASE}"

  # Reads are single-end, compressed, and aligned to a pre-built index.  The
  # transcriptome BAM and GeneCounts outputs are needed by downstream steps.
  STAR \
    --runThreadN '"${THREADS_PER_JOB}"' \
    --genomeDir "'"${STAR_INDEX}"'" \
    --readFilesIn "$FQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${PREFIX}" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --quantTranscriptomeSAMoutput BanSingleEnd \
    --twopassMode Basic \
    --alignIntronMax 2500 \
    --limitBAMsortRAM 100000000000 \
    > "'"${OUTDIR}"'/${BASE}.STAR.log" 2>&1

  # Indexing is convenient for downstream random access.  Keep it optional so
  # STAR alignment itself can still complete on systems without samtools.
  if command -v samtools >/dev/null 2>&1; then
    samtools index -@ 4 "${PREFIX}Aligned.sortedByCoord.out.bam"
  fi
' _ {}

echo "Done.
- Output dir:        ${OUTDIR}
- Per-sample logs:   ${OUTDIR}/*.STAR.log
- Coord BAM:         ${OUTDIR}/*Aligned.sortedByCoord.out.bam
- Transcriptome BAM: ${OUTDIR}/*Aligned.toTranscriptome.out.bam
- GeneCounts:        ${OUTDIR}/*ReadsPerGene.out.tab"
