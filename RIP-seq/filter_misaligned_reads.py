#!/usr/bin/env python3
"""Identify misaligned reads from STAR transcriptome SAM records.

This script intentionally avoids non-standard Python dependencies. Feed it
`samtools view` output on stdin. It writes read names to remove and a summary
CSV. The shell wrapper uses those read names to create the clean/discarded BAMs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


# Coordinates in this script are zero-based, half-open intervals.  Keeping the
# annotation records explicit makes the projection logic readable and avoids
# requiring non-standard Python packages.
@dataclass
class Gene:
    gene_id: str
    gff_id: str
    seqid: str
    start: int
    end: int
    strand: str
    biotype: str = ""


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    segments: list[dict] = field(default_factory=list)


def parse_attrs(attr_text: str) -> dict[str, str]:
    """Parse GFF3 key=value and GTF-like key value attribute fields."""
    attrs: dict[str, str] = {}
    for field in attr_text.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value.strip().strip('"')
        elif " " in field:
            key, value = field.split(" ", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def parse_gff(path: Path) -> tuple[dict[str, Gene], dict[str, Transcript]]:
    """Read gene/transcript/exon models and build transcript-coordinate maps."""
    genes_by_gff_id: dict[str, Gene] = {}
    genes_by_locus: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}

    # GFF headers and malformed records are ignored; feature records are
    # linked in file order because children refer to previously seen parents.
    with path.open(newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature, start_s, end_s, _score, strand, _phase, attrs_s = fields
            attrs = parse_attrs(attrs_s)
            start = int(start_s) - 1
            end = int(end_s)

            if feature == "gene":
                # locus_tag is the stable identifier used by the surrounding
                # workflow; the GFF ID is the fallback when it is absent.
                gff_id = attrs.get("ID", "")
                gene_id = attrs.get("locus_tag") or gff_id.removeprefix("gene-")
                if not gff_id or not gene_id:
                    continue
                gene = Gene(
                    gene_id=gene_id,
                    gff_id=gff_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    biotype=attrs.get("gene_biotype", ""),
                )
                genes_by_gff_id[gff_id] = gene
                genes_by_locus[gene_id] = gene
                continue

            if feature in {"mRNA", "lnc_RNA", "rRNA", "snRNA", "snoRNA", "tRNA", "SRP_RNA"}:
                # Only transcripts with a known parent gene participate in the
                # later overlap tests.
                tid = attrs.get("ID", "")
                parent = attrs.get("Parent", "")
                parent_gene = genes_by_gff_id.get(parent)
                if tid and parent_gene:
                    transcripts[tid] = Transcript(
                        transcript_id=tid,
                        gene_id=parent_gene.gene_id,
                        seqid=seqid,
                        strand=strand,
                    )
                continue

            if feature == "exon":
                # Exons are collected first and converted to transcript
                # segments after all annotation lines have been read.
                parent = attrs.get("Parent", "")
                tx = transcripts.get(parent)
                if tx:
                    tx.exons.append((start, end))

    # Build cumulative transcript coordinates. Reverse-strand exons are sorted
    # in transcript direction before their genomic mapping is recorded.
    for tx in transcripts.values():
        ordered = sorted(tx.exons, key=lambda item: item[0], reverse=(tx.strand == "-"))
        tx_pos = 0
        for g_start, g_end in ordered:
            length = g_end - g_start
            tx.segments.append(
                {
                    "tx_start": tx_pos,
                    "tx_end": tx_pos + length,
                    "g_start": g_start,
                    "g_end": g_end,
                }
            )
            tx_pos += length

    return genes_by_locus, transcripts


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Return whether two half-open intervals overlap."""
    return a_start < b_end and b_start < a_end


def interval_contains(outer_start: int, outer_end: int, inner_start: int, inner_end: int) -> bool:
    """Return whether the inner interval is fully contained by the outer one."""
    return outer_start <= inner_start and inner_end <= outer_end


def blocks_within_gene_with_tolerance(
    blocks: list[tuple[int, int]], gene: Gene, tolerance_nt: int
) -> bool:
    """Check that every aligned block lies within a gene plus edge tolerance."""
    return bool(blocks) and all(
        interval_contains(
            gene.start - tolerance_nt,
            gene.end + tolerance_nt,
            block_start,
            block_end,
        )
        for block_start, block_end in blocks
    )


def parse_cigar_blocks(pos_1based: int, cigar: str) -> list[tuple[int, int]]:
    """Extract reference-consuming alignment blocks from a SAM CIGAR string."""
    pos = pos_1based - 1
    blocks: list[tuple[int, int]] = []
    for length_s, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_s)
        if op in {"M", "=", "X"}:
            blocks.append((pos, pos + length))
            pos += length
        elif op in {"D", "N"}:
            pos += length
        elif op in {"I", "S", "H", "P"}:
            continue
    return blocks


def map_tx_blocks_to_genome(
    tx: Transcript, tx_blocks: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Project transcript-coordinate blocks onto the transcript's genomic exons."""
    genomic_blocks: list[tuple[int, int]] = []
    for block_start, block_end in tx_blocks:
        for seg in tx.segments:
            overlap_start = max(block_start, seg["tx_start"])
            overlap_end = min(block_end, seg["tx_end"])
            if overlap_start >= overlap_end:
                continue
            rel_start = overlap_start - seg["tx_start"]
            rel_end = overlap_end - seg["tx_start"]
            if tx.strand == "-":
                g_start = seg["g_end"] - rel_end
                g_end = seg["g_end"] - rel_start
            else:
                g_start = seg["g_start"] + rel_start
                g_end = seg["g_start"] + rel_end
            genomic_blocks.append((g_start, g_end))
    return genomic_blocks


def collect_same_strand_overlapping_unintended_genes(
    full_genes: dict[str, Gene], reference_genes: dict[str, Gene]
) -> dict[str, list[Gene]]:
    """Index unintended genes overlapping each reference gene on the same strand."""
    overlapping_by_reference_gene: dict[str, list[Gene]] = defaultdict(list)
    unintended = [gene for gid, gene in full_genes.items() if gid not in reference_genes]
    for ref_gene in reference_genes.values():
        for unintended_gene in unintended:
            if unintended_gene.seqid != ref_gene.seqid:
                continue
            if unintended_gene.strand != ref_gene.strand:
                continue
            if interval_overlap(ref_gene.start, ref_gene.end, unintended_gene.start, unintended_gene.end):
                overlapping_by_reference_gene[ref_gene.gene_id].append(unintended_gene)
    return overlapping_by_reference_gene


def parse_locus_list(loci_text: str) -> set[str]:
    """Convert a comma-separated locus argument into a trimmed set."""
    return {item.strip() for item in loci_text.split(",") if item.strip()}


def main() -> int:
    # The shell wrapper supplies SAM on stdin and names the two output files;
    # these options expose the annotation and filtering policy for reuse.
    parser = argparse.ArgumentParser(description="Identify transcriptome-BAM reads to remove.")
    parser.add_argument("--full-gff", required=True, type=Path)
    parser.add_argument("--reference-gff", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--removed-reads", required=True, type=Path)
    parser.add_argument("--summary-csv", required=True, type=Path)
    parser.add_argument("--special-locus", default="cpbgf_300675")
    parser.add_argument("--always-remove-loci", default="cpbgf_3001410")
    parser.add_argument("--edge-tolerance-nt", default=5, type=int)
    parser.add_argument(
        "--same-orientation-read-strand",
        choices=("reverse", "forward"),
        default="reverse",
        help="Expected read orientation relative to the transcript for the same-orientation rule (SR=reverse; SF=forward).",
    )
    args = parser.parse_args()
    if args.edge_tolerance_nt < 0:
        raise SystemExit("ERROR: --edge-tolerance-nt must be >= 0")
    always_remove_loci = parse_locus_list(args.always_remove_loci)

    # The full annotation identifies unintended overlaps. The restricted GFF
    # defines the transcript models that STAR was intended to align against.
    full_genes, _full_transcripts = parse_gff(args.full_gff)
    reference_genes, reference_transcripts = parse_gff(args.reference_gff)
    overlapping_by_reference_gene = collect_same_strand_overlapping_unintended_genes(
        full_genes, reference_genes
    )

    # Fail early when configured loci do not exist in the selected references;
    # silently accepting a typo would produce an apparently valid but wrong BAM.
    special_gene = full_genes.get(args.special_locus)
    if special_gene is None:
        raise SystemExit(f"ERROR: special locus not found in full GFF: {args.special_locus}")
    missing_always_remove_loci = sorted(
        locus for locus in always_remove_loci if locus not in reference_genes
    )
    if missing_always_remove_loci:
        raise SystemExit(
            "ERROR: always-remove locus not found in reference GFF: "
            + ",".join(missing_always_remove_loci)
        )

    # These are read-name sets. A read with multiple transcript records is
    # removed consistently from both BAMs once its name enters any set.
    always_remove_reads: set[str] = set()
    special_reads: set[str] = set()
    same_strand_overlap_reads: set[str] = set()
    same_strand_overlap_trigger_genes_by_read: dict[str, set[str]] = defaultdict(set)
    mapped_genes_by_read: dict[str, set[str]] = defaultdict(set)
    unmodeled_reference_names: set[str] = set()
    reference_names_without_exons: set[str] = set()
    stats = defaultdict(int)

    # Stream SAM records one at a time so memory usage does not scale with BAM
    # size. Header lines are ignored because they carry no read evidence.
    for line in sys.stdin:
        if not line.strip() or line.startswith("@"):
            continue
        stats["tx_records_seen"] += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            stats["malformed_sam_records"] += 1
            continue
        qname = fields[0]
        flag = int(fields[1])
        refname = fields[2]
        pos = int(fields[3])
        cigar = fields[5]
        # Unmapped records cannot be projected to a gene and therefore cannot
        # trigger a cleanup rule.
        if flag & 4 or refname == "*" or cigar == "*":
            continue
        # Every mapped transcript name must have a model in the reference GFF;
        # unknown names are collected for a clear error after streaming ends.
        tx = reference_transcripts.get(refname)
        if tx is None:
            stats["tx_records_without_reference_transcript_model"] += 1
            unmodeled_reference_names.add(refname)
            continue
        if not tx.segments:
            stats["tx_records_without_exon_model"] += 1
            reference_names_without_exons.add(refname)
            continue
        stats["tx_records_inspected"] += 1
        mapped_genes_by_read[qname].add(tx.gene_id)
        if tx.gene_id in always_remove_loci:
            always_remove_reads.add(qname)

        # Convert the transcriptome CIGAR to transcript blocks, then project
        # those blocks back across the transcript's exon structure.
        tx_blocks = parse_cigar_blocks(pos, cigar)
        genomic_blocks = map_tx_blocks_to_genome(tx, tx_blocks)
        if not genomic_blocks:
            continue

        # Direct special-locus overlap is deliberately evaluated before strand
        # filtering because it is an unconditional cleanup rule.
        if tx.seqid == special_gene.seqid and any(
            interval_overlap(start, end, special_gene.start, special_gene.end)
            for start, end in genomic_blocks
        ):
            special_reads.add(qname)

        # Apply stranded-library compatibility only to the same-orientation
        # cleanup. The direct cpbgf_3001410 and cpbgf_300675 rules above
        # intentionally inspect every mapped record.
        is_reverse = bool(flag & 16)
        if args.same_orientation_read_strand == "reverse" and not is_reverse:
            continue
        if args.same_orientation_read_strand == "forward" and is_reverse:
            continue

        # A projected read fully inside an unintended same-strand gene is
        # ambiguous and is marked for removal; retain trigger IDs for auditing.
        for unintended_gene in overlapping_by_reference_gene.get(tx.gene_id, []):
            if blocks_within_gene_with_tolerance(
                genomic_blocks, unintended_gene, args.edge_tolerance_nt
            ):
                same_strand_overlap_reads.add(qname)
                same_strand_overlap_trigger_genes_by_read[qname].add(unintended_gene.gene_id)

    if stats["malformed_sam_records"]:
        raise SystemExit(
            f"ERROR: {stats['malformed_sam_records']} malformed SAM records were encountered."
        )
    if unmodeled_reference_names:
        preview = ",".join(sorted(unmodeled_reference_names)[:10])
        suffix = "" if len(unmodeled_reference_names) <= 10 else ",..."
        raise SystemExit(
            "ERROR: transcriptome BAM reference names were not found in the reference GFF "
            f"({len(unmodeled_reference_names)} unique; examples: {preview}{suffix})."
        )
    if reference_names_without_exons:
        preview = ",".join(sorted(reference_names_without_exons)[:10])
        suffix = "" if len(reference_names_without_exons) <= 10 else ",..."
        raise SystemExit(
            "ERROR: reference-GFF transcripts used by the BAM have no exon model "
            f"({len(reference_names_without_exons)} unique; examples: {preview}{suffix})."
        )

    # Make the summary mutually exclusive by assigning each read to the first
    # applicable cleanup step, while retaining the union as the final removal set.
    always_remove_step_reads = always_remove_reads
    special_step_reads = special_reads - always_remove_step_reads
    same_strand_overlap_step_reads = same_strand_overlap_reads - always_remove_step_reads - special_step_reads
    all_removed_reads = always_remove_reads | special_reads | same_strand_overlap_reads
    same_strand_overlap_step_trigger_genes = {
        gene_id
        for read_name in same_strand_overlap_step_reads
        for gene_id in same_strand_overlap_trigger_genes_by_read.get(read_name, set())
    }
    all_same_strand_overlap_trigger_genes = {
        gene_id
        for read_name in same_strand_overlap_reads
        for gene_id in same_strand_overlap_trigger_genes_by_read.get(read_name, set())
    }

    # Write sorted names for deterministic output and reproducible downstream
    # filtering. Parent directories are created only after validation succeeds.
    args.removed_reads.parent.mkdir(parents=True, exist_ok=True)
    with args.removed_reads.open("w") as handle:
        for read_name in sorted(all_removed_reads):
            handle.write(f"{read_name}\n")

    def mapped_gene_list(reads: set[str]) -> str:
        genes: set[str] = set()
        for read_name in reads:
            genes.update(mapped_genes_by_read.get(read_name, set()))
        return ";".join(sorted(genes))

    def gene_list(genes: set[str]) -> str:
        return ";".join(sorted(genes))

    # One row per cleanup stage records counts, trigger loci, and mapped genes.
    rows = [
        {
            "sample": args.sample,
            "cleanup_step": "always_remove_aligned_loci",
            "reads_removed": len(always_remove_step_reads),
            "trigger_gene_ids": gene_list(always_remove_loci if always_remove_step_reads else set()),
            "mapped_gene_ids_before_removal": mapped_gene_list(always_remove_step_reads),
        },
        {
            "sample": args.sample,
            "cleanup_step": f"{args.special_locus}_overlap",
            "reads_removed": len(special_step_reads),
            "trigger_gene_ids": args.special_locus if special_step_reads else "",
            "mapped_gene_ids_before_removal": mapped_gene_list(special_step_reads),
        },
        {
            "sample": args.sample,
            "cleanup_step": "same_orientation_overlapping_unintended",
            "reads_removed": len(same_strand_overlap_step_reads),
            "trigger_gene_ids": gene_list(same_strand_overlap_step_trigger_genes),
            "mapped_gene_ids_before_removal": mapped_gene_list(same_strand_overlap_step_reads),
        },
        {
            "sample": args.sample,
            "cleanup_step": "total_unique_removed",
            "reads_removed": len(all_removed_reads),
            "trigger_gene_ids": gene_list(
                (always_remove_loci if always_remove_step_reads else set())
                |
                ({args.special_locus} if special_reads else set())
                | all_same_strand_overlap_trigger_genes
            ),
            "mapped_gene_ids_before_removal": mapped_gene_list(all_removed_reads),
        },
    ]

    # Repeat the run metadata on every row so the CSV remains self-describing.
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.summary_csv.open("w", newline="") as handle:
        fieldnames = [
            "sample",
            "cleanup_step",
            "reads_removed",
            "trigger_gene_ids",
            "mapped_gene_ids_before_removal",
            "tx_records_seen",
            "tx_records_inspected",
            "tx_records_without_reference_transcript_model",
            "tx_records_without_exon_model",
            "malformed_sam_records",
            "edge_tolerance_nt",
            "always_remove_loci",
            "same_orientation_read_strand",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "tx_records_seen": stats["tx_records_seen"],
                    "tx_records_inspected": stats["tx_records_inspected"],
                    "tx_records_without_reference_transcript_model": stats[
                        "tx_records_without_reference_transcript_model"
                    ],
                    "tx_records_without_exon_model": stats[
                        "tx_records_without_exon_model"
                    ],
                    "malformed_sam_records": stats["malformed_sam_records"],
                    "edge_tolerance_nt": args.edge_tolerance_nt,
                    "always_remove_loci": gene_list(always_remove_loci),
                    "same_orientation_read_strand": args.same_orientation_read_strand,
                }
            )

    print(f"Removed reads: {len(all_removed_reads)}", file=sys.stderr)
    print(f"Wrote: {args.removed_reads}", file=sys.stderr)
    print(f"Wrote: {args.summary_csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
