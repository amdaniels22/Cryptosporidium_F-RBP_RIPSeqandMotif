#!/usr/bin/env bash
# Identify and remove transcriptome alignments whose genomic projection is
# likely confounded by unintended overlapping loci.  The Python program reads
# SAM records on stdin; this wrapper handles BAM discovery, parallelism, logs,
# and construction of cleaned/discarded BAM files.
set -euo pipefail

# Clean STAR transcriptome BAMs before downstream quantification.
#
# Default run location:
#   Run this from the RIP-seq project directory. The script looks for:
#     ./star_out_mRNA_lncRNA_only_clean_500bp/*.Aligned.toTranscriptome.out.bam
#     ./star_out_mRNA_lncRNA_only_clean_500bp/*.Aligned.sortedByCoord.out.bam

# These environment-overridable settings control parallel jobs and samtools
# indexing threads without requiring edits to the script.
JOBS="${JOBS:-4}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"

# Folder containing matching STAR BAM pairs.  For every transcriptome BAM the
# wrapper expects a coordinate-sorted BAM with the same sample prefix.
#   *.Aligned.toTranscriptome.out.bam
#   *.Aligned.sortedByCoord.out.bam
BAM_DIR="${BAM_DIR:-${STAR_DIR:-./star_out_mRNA_lncRNA_only_clean_500bp}}"
OUTDIR="${OUTDIR:-./misaligned_read_filter_500bp}"
LOG_DIR="${LOG_DIR:-${OUTDIR}/logs}"

FULL_GFF="${FULL_GFF:-./CpBGF_genome_V1.1.gff}"
REFERENCE_GFF="${REFERENCE_GFF:-./CpBGF_genome_V1.1_mRNA_lncRNA_only_clean.gff}"
SPECIAL_LOCUS="${SPECIAL_LOCUS:-cpbgf_300675}"
ALWAYS_REMOVE_LOCI="${ALWAYS_REMOVE_LOCI:-cpbgf_3001410}"
EDGE_TOLERANCE_NT="${EDGE_TOLERANCE_NT:-5}"
# Salmon library type SR: valid reads are reverse relative to the transcript.
# Use SAME_ORIENTATION_READ_STRAND=forward for an SF library.
SAME_ORIENTATION_READ_STRAND="${SAME_ORIENTATION_READ_STRAND:-reverse}"

# Resolve the Python script relative to this wrapper, so it works regardless of
# the caller's current directory (the input/output defaults remain project-
# relative as documented in the README).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER_PY="${FILTER_PY:-${SCRIPT_DIR}/filter_misaligned_reads.py}"

mkdir -p "${OUTDIR}" "${LOG_DIR}"

# Fail early with a useful message instead of starting workers that will all
# fail for the same missing dependency or reference file.
for required in samtools "${PYTHON}" "${FILTER_PY}" "${FULL_GFF}" "${REFERENCE_GFF}"; do
  if [[ "${required}" == "samtools" ]]; then
    command -v samtools >/dev/null 2>&1 || {
      echo "ERROR: samtools not found on PATH." >&2
      exit 1
    }
  elif [[ "${required}" == "${PYTHON}" ]]; then
    command -v "${PYTHON}" >/dev/null 2>&1 || {
      echo "ERROR: ${PYTHON} not found on PATH." >&2
      exit 1
    }
  elif [[ ! -e "${required}" ]]; then
    echo "ERROR: required file not found: ${required}" >&2
    exit 1
  fi
done

# Export values so the nested xargs shell receives the same configuration.
export ALWAYS_REMOVE_LOCI EDGE_TOLERANCE_NT FILTER_PY FULL_GFF LOG_DIR OUTDIR PYTHON REFERENCE_GFF SAME_ORIENTATION_READ_STRAND SAMTOOLS_THREADS SPECIAL_LOCUS

shopt -s nullglob
# nullglob makes an unmatched glob expand to no elements rather than to a
# literal filename.
TX_BAMS=( "${BAM_DIR}"/*.Aligned.toTranscriptome.out.bam )
shopt -u nullglob

if [[ ${#TX_BAMS[@]} -eq 0 ]]; then
  echo "ERROR: no *.Aligned.toTranscriptome.out.bam files found in ${BAM_DIR}." >&2
  exit 1
fi

# NUL-delimit filenames so spaces or other shell metacharacters are safe.
printf '%s\0' "${TX_BAMS[@]}" \
| xargs -0 -n1 -P "${JOBS}" bash -c '
  set -euo pipefail

  # $1 is the filename supplied by xargs; derive all related output names from
  # the common STAR sample prefix.
  tx_bam="$1"
  base="${tx_bam%.Aligned.toTranscriptome.out.bam}"
  sample="$(basename "${base}")"
  coord_bam="${base}.Aligned.sortedByCoord.out.bam"
  out_prefix="${OUTDIR}/${sample}"
  clean_tx="${out_prefix}.Aligned.toTranscriptome.clean.out.bam"
  discarded_coord="${out_prefix}.Aligned.sortedByCoord.discarded.out.bam"
  removed_reads="${out_prefix}.removed_read_names.txt"
  summary_csv="${out_prefix}.cleanup_summary.csv"

  if [[ ! -f "${coord_bam}" ]]; then
    echo "ERROR: missing coordinate BAM for ${sample}: ${coord_bam}" >&2
    exit 1
  fi

  echo "[clean] ${sample}"
  # Stream the transcriptome BAM through Python.  stdout is the diagnostic log
  # and stderr is captured separately for troubleshooting.
  samtools view "${tx_bam}" \
    | "${PYTHON}" "${FILTER_PY}" \
        --full-gff "${FULL_GFF}" \
        --reference-gff "${REFERENCE_GFF}" \
        --special-locus "${SPECIAL_LOCUS}" \
        --always-remove-loci "${ALWAYS_REMOVE_LOCI}" \
        --edge-tolerance-nt "${EDGE_TOLERANCE_NT}" \
        --same-orientation-read-strand "${SAME_ORIENTATION_READ_STRAND}" \
        --sample "${sample}" \
        --removed-reads "${removed_reads}" \
        --summary-csv "${summary_csv}" \
        > "${LOG_DIR}/${sample}.filter.stdout.log" \
        2> "${LOG_DIR}/${sample}.filter.stderr.log"

  # If any read names were flagged, remove them from the transcriptome BAM and
  # retain exactly those reads in a separate coordinate BAM for auditing.
  if [[ -s "${removed_reads}" ]]; then
    samtools view -h "${tx_bam}" \
      | awk "NR==FNR {drop[\$1]=1; next} /^@/ {print; next} !(\$1 in drop)" "${removed_reads}" - \
      | samtools view -b -o "${clean_tx}"

    samtools view -h "${coord_bam}" \
      | awk "NR==FNR {drop[\$1]=1; next} /^@/ {print; next} (\$1 in drop)" "${removed_reads}" - \
      | samtools view -b -o "${discarded_coord}"
  else
    # Preserve the header even when no alignments are discarded, creating a
    # valid empty BAM rather than an absent output.
    echo "[clean] ${sample}: no reads met the removal criteria"
    samtools view -h "${tx_bam}" \
      | samtools view -b -o "${clean_tx}" -
    samtools view -H "${coord_bam}" \
      | samtools view -b -o "${discarded_coord}" -
  fi

  # quickcheck catches truncated BAMs before they propagate downstream; only
  # the discarded coordinate BAM is indexed because it is coordinate-sorted.
  samtools quickcheck "${clean_tx}" "${discarded_coord}"
  samtools index -@ "${SAMTOOLS_THREADS}" "${discarded_coord}"
' _

echo "Done.
- BAM input dir:             ${BAM_DIR}
- Output dir:                ${OUTDIR}
- Clean transcriptome BAMs:  ${OUTDIR}/*.Aligned.toTranscriptome.clean.out.bam
- Discarded coordinate BAMs: ${OUTDIR}/*.Aligned.sortedByCoord.discarded.out.bam
- Discarded BAM indexes:     ${OUTDIR}/*.Aligned.sortedByCoord.discarded.out.bam.bai
- Summary CSVs:              ${OUTDIR}/*.cleanup_summary.csv
- Removed read lists:        ${OUTDIR}/*.removed_read_names.txt
- Logs:                      ${LOG_DIR}"
