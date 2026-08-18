#!/usr/bin/env bash
# Project the final set of retained, genomically unique read names back onto
# coordinate-sorted STAR BAMs.  These BAMs are used to make genome-browser and
# coverage tracks while preserving the same read set used for Salmon.
set -euo pipefail
shopt -s nullglob

# -------- Settings --------
COORD_DIR="${COORD_DIR:-./star_out_mRNA_lncRNA_only_clean_500bp}"
UNIQUE_TX_DIR="${UNIQUE_TX_DIR:-./unique_transcriptome_bams_500bp}"
OUTDIR="${OUTDIR:-./coord_bams_500_clean_unique}"
INDEX_THREADS=4

mkdir -p "$OUTDIR"

# -------- Requirements --------
for program in samtools cut sort; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

# -------- Locate coordinate BAMs --------
# With nullglob enabled, a missing match yields an empty array and a clear
# error below instead of accidentally processing a literal glob.
coord_bams=("$COORD_DIR"/*Aligned.sortedByCoord.out.bam)

if (( ${#coord_bams[@]} == 0 )); then
    echo "ERROR: No coordinate BAMs found in: $COORD_DIR" >&2
    exit 1
fi

# -------- Process each sample --------
for coord_bam in "${coord_bams[@]}"; do
    # Derive the sample name from the STAR coordinate-BAM suffix and use it to
    # locate the matching unique transcriptome BAM.
    filename=$(basename "$coord_bam")
    base=${filename%.Aligned.sortedByCoord.out.bam}

    transcript_bam="$UNIQUE_TX_DIR/${base}.Aligned.toTranscriptome.clean.unique.out.bam"
    read_names="$OUTDIR/${base}.retained_transcriptome.read_names.txt"
    output_bam="$OUTDIR/${base}.Aligned.sortedByCoord.clean.unique.out.bam"

    echo "[Sample] $base"

    if [[ ! -f "$transcript_bam" ]]; then
        echo "ERROR: Missing unique transcriptome BAM:" >&2
        echo "       $transcript_bam" >&2
        exit 1
    fi

    if ! samtools quickcheck "$coord_bam" "$transcript_bam"; then
        echo "ERROR: Invalid or truncated input BAM for $base" >&2
        exit 1
    fi

    # The unique transcriptome BAM is the definitive retained-read set. Do not
    # inspect transcriptome NH here: it can reflect alternate transcript
    # projections rather than genomic multimapping.
    samtools view "$transcript_bam" \
        | cut -f1 \
        | LC_ALL=C sort -u \
        > "$read_names"

    retained_reads=$(wc -l < "$read_names")
    retained_reads=${retained_reads//[[:space:]]/}

    if [[ "$retained_reads" -eq 0 ]]; then
        echo "WARNING: No retained transcriptome reads found for $base" >&2

        # Create a valid empty BAM with the original coordinate-BAM header.
        samtools view -H "$coord_bam" \
            | samtools view -b -o "$output_bam" -
    else
        # Filter the coordinate BAM by the retained names, NH:i:1, and primary
        # mapped flags.  This excludes unmapped, secondary, and supplementary
        # records while retaining the exact final genomic read set.
        # Keep coordinate alignments whose read names occur in the final unique
        # transcriptome BAM and that are:
        #   - uniquely mapped: NH:i:1
        #   - mapped and primary: exclude flags 4, 256 and 2048
        samtools view \
            -b \
            -N "$read_names" \
            -d NH:1 \
            -F 2308 \
            -o "$output_bam" \
            "$coord_bam"
    fi

    samtools index -@ "$INDEX_THREADS" "$output_bam"
    samtools quickcheck -v "$output_bam"

    output_records=$(samtools view -c "$output_bam")
    output_reads=$(samtools view "$output_bam" | cut -f1 | LC_ALL=C sort -u | wc -l)
    output_reads=${output_reads//[[:space:]]/}

    if [[ "$output_reads" -ne "$retained_reads" ]]; then
        echo "ERROR: Coordinate/transcriptome retained-read mismatch for $base " >&2
        echo "       transcriptome=$retained_reads coordinate=$output_reads" >&2
        exit 1
    fi

    echo "  Retained transcriptome read names: $retained_reads"
    echo "  Retained coordinate read names:    $output_reads"
    echo "  Coordinate BAM records:      $output_records"
    echo "  Output: $output_bam"
done

echo
echo "Done."
echo "Output directory: $OUTDIR"
