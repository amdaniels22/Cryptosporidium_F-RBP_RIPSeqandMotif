#!/usr/bin/env bash
# Keep transcriptome alignments belonging to reads that are uniquely mapped in
# the genomic (coordinate-sorted) STAR BAM.  This separates genomic uniqueness
# from transcript-level projections, which may legitimately have alternatives.
set -euo pipefail
shopt -s nullglob

# Clean transcriptome BAMs produced by filter_misaligned_reads.sh.
# Override these paths with INPUT_DIR, COORD_DIR, or OUTPUT_DIR when running a
# differently named pipeline.
INPUT_DIR="${INPUT_DIR:-./misaligned_read_filter_500bp}"
# Original coordinate BAMs used to determine genomic uniqueness.
COORD_DIR="${COORD_DIR:-./star_out_mRNA_lncRNA_only_clean_500bp}"
OUTPUT_DIR="${OUTPUT_DIR:-./unique_transcriptome_bams_500bp}"

mkdir -p "$OUTPUT_DIR"

# nullglob ensures an empty input directory is detected by the explicit check.
bams=("$INPUT_DIR"/*Aligned.toTranscriptome.clean.out.bam)

if (( ${#bams[@]} == 0 )); then
    echo "ERROR: No cleaned transcriptome BAMs found in $INPUT_DIR" >&2
    exit 1
fi

# samtools performs BAM operations; cut/sort construct a unique read-name list.
for program in samtools cut sort; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

for bam in "${bams[@]}"; do
    # All output names retain the STAR sample prefix, making them match the
    # coordinate-BAM stage that follows.
    filename=$(basename "$bam")
    base=${filename%.Aligned.toTranscriptome.clean.out.bam}
    coord_bam="$COORD_DIR/${base}.Aligned.sortedByCoord.out.bam"
    read_names="$OUTPUT_DIR/${base}.coordinate_unique.read_names.txt"
    output="$OUTPUT_DIR/${base}.Aligned.toTranscriptome.clean.unique.out.bam"

    echo "Filtering $base"

    if [[ ! -f "$coord_bam" ]]; then
        echo "ERROR: Missing coordinate BAM: $coord_bam" >&2
        exit 1
    fi

    if ! samtools quickcheck "$bam" "$coord_bam"; then
        echo "ERROR: Invalid or truncated input BAM for $base" >&2
        exit 1
    fi

    # STAR NH in the transcriptome BAM can count alternate transcript
    # projections. Determine uniqueness only from mapped primary genomic records.
    samtools view -F 2308 -d NH:1 "$coord_bam" \
        | cut -f1 \
        | LC_ALL=C sort -u \
        > "$read_names"

    # -N selects reads by name while preserving every transcriptome record for
    # each selected read, including alternate transcript projections.
    if [[ -s "$read_names" ]]; then
        # Preserve every transcriptome record for each genomically unique read,
        # including alternate transcript projections needed by Salmon.
        samtools view -b -N "$read_names" -o "$output" "$bam"
    else
        # A header-only BAM is still a valid BAM and makes the empty result
        # explicit for downstream tools.
        samtools view -H "$bam" | samtools view -b -o "$output" -
    fi

    # Validate the output before reporting counts.
    samtools quickcheck -v "$output"

    unique_reads=$(wc -l < "$read_names")
    unique_reads=${unique_reads//[[:space:]]/}
    retained_reads=$(samtools view "$output" | cut -f1 | LC_ALL=C sort -u | wc -l)
    retained_reads=${retained_reads//[[:space:]]/}
    echo "  Input records:  $(samtools view -c "$bam")"
    echo "  Genomically unique coordinate read names: $unique_reads"
    echo "  Retained transcriptome read names: $retained_reads"
    echo "  Output transcriptome records: $(samtools view -c "$output")"
done
