#!/usr/bin/env python3
"""
Filter overlapping genes with priority for protein-coding genes.
Strategy:
1. If overlap includes protein-coding gene(s), keep only protein-coding
2. Among genes of same type, keep the one with most exons
3. Generate lists of genes to keep and remove

Usage: python filter_overlapping_genes.py <input.gff> <gtf_directory> [output.gff]
"""

import sys
import os
import collections
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Filter overlapping genes with priority for protein-coding genes.',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Filtering strategy:
  1. If overlap includes protein-coding gene(s), keep only protein-coding
  2. Among genes of same type, keep the one with most exons
  3. Generate two text files listing genes to keep and remove

Example:
  python filter_overlapping_genes.py gencode.v28.dexseq.bygene.cat.gff gtf/ overlapping_genes
"""
)
parser.add_argument('input_gff', help='Input GFF file with all genes')
parser.add_argument('gtf_dir', help='Directory containing per-gene GTF files (gene_id.gtf)')
parser.add_argument('output_prefix', nargs='?', default='overlapping_genes',
                    help='Prefix for output files (default: overlapping_genes)')

args = parser.parse_args()

# Set up file paths
input_gff = args.input_gff
gtf_dir = args.gtf_dir
output_prefix = args.output_prefix

# Validate inputs
if not os.path.exists(input_gff):
    print(f"Error: Input GFF file not found: {input_gff}", file=sys.stderr)
    sys.exit(1)

if not os.path.isdir(gtf_dir):
    print(f"Error: GTF directory not found: {gtf_dir}", file=sys.stderr)
    sys.exit(1)

print(f"Input GFF: {input_gff}", file=sys.stderr)
print(f"GTF directory: {gtf_dir}", file=sys.stderr)
print(f"Output prefix: {output_prefix}", file=sys.stderr)
print("", file=sys.stderr)

def get_gene_biotype(gene_id, gtf_directory):
    """Get gene biotype from GTF file"""
    gtf_file = os.path.join(gtf_directory, f"{gene_id}.gtf")
    try:
        with open(gtf_file) as f:
            for line in f:
                if line.startswith('chr') and '\tgene\t' in line:
                    if 'gene_type "' in line:
                        biotype = line.split('gene_type "')[1].split('"')[0]
                        return biotype
    except:
        return "unknown"
    return "unknown"

print("Reading gene coordinates and exonic parts...", file=sys.stderr)
gene_coords = {}  # gene_id -> (chrom, start, end, strand)
exon_counts = {}  # gene_id -> number of exonic_parts
exonic_parts = []  # list of (chrom, start, end, strand, gene_id)

with open(input_gff) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        
        if 'gene_id "' not in parts[8]:
            continue
        gene_id = parts[8].split('gene_id "')[1].split('"')[0]
        
        if parts[2] == 'aggregate_gene':
            gene_coords[gene_id] = (parts[0], int(parts[3]), int(parts[4]), parts[6])
            exon_counts[gene_id] = 0
        elif parts[2] == 'exonic_part':
            exon_counts[gene_id] = exon_counts.get(gene_id, 0) + 1
            # Store exonic part info for overlap detection
            exonic_parts.append((parts[0], int(parts[3]), int(parts[4]), parts[6], gene_id))

print(f"Found {len(gene_coords)} genes with {len(exonic_parts)} exonic parts", file=sys.stderr)

# Read gene biotypes
print("Reading gene biotypes from GTF files...", file=sys.stderr)
gene_biotypes = {}
for i, gene_id in enumerate(gene_coords.keys()):
    if i % 5000 == 0:
        print(f"  Processed {i}/{len(gene_coords)} genes...", file=sys.stderr)
    gene_biotypes[gene_id] = get_gene_biotype(gene_id, gtf_dir)

print(f"Retrieved biotypes for {len(gene_biotypes)} genes", file=sys.stderr)

# Count biotypes
biotype_counts = collections.Counter(gene_biotypes.values())
print("\nBiotype distribution:", file=sys.stderr)
for biotype, count in biotype_counts.most_common(10):
    print(f"  {biotype}: {count}", file=sys.stderr)

# Detect overlapping genes using same logic as dexseq_prepare_annotation.py
# Two genes overlap if they share exonic bases on the same strand
print("\nDetecting overlaps (exon-level, stranded)...", file=sys.stderr)

# Group exonic parts by chromosome and strand
chrom_strand_exons = collections.defaultdict(list)
for chrom, start, end, strand, gene_id in exonic_parts:
    chrom_strand_exons[(chrom, strand)].append((start, end, gene_id))

# For each chromosome+strand, find overlapping exons and their genes
# This mimics HTSeq.GenomicArrayOfSets behavior
gene_sets = collections.defaultdict(lambda: set())  # gene_id -> set of genes it overlaps with

for (chrom, strand), exons in chrom_strand_exons.items():
    # Sort by start position
    exons.sort(key=lambda x: x[0])
    
    # Create intervals and track which genes are in each
    # Use sweep line algorithm to find overlapping exons
    events = []  # (position, type, gene_id) where type is 'start' or 'end'
    for start, end, gene_id in exons:
        events.append((start, 'start', gene_id))
        events.append((end + 1, 'end', gene_id))  # +1 because end is inclusive
    
    events.sort()
    
    # Track active genes at each position
    active_genes = set()
    for pos, event_type, gene_id in events:
        if event_type == 'start':
            # Gene starts: it overlaps with all currently active genes
            if active_genes:
                # Genes that share exonic bases
                full_set = set([gene_id])
                for active_gene in active_genes:
                    full_set.add(active_gene)
                    # Also add genes that previously overlapped with active genes
                    full_set |= gene_sets[active_gene]
                
                # Update all genes in full_set to know about each other
                for g in full_set:
                    gene_sets[g] = full_set
            else:
                # First gene at this position
                if gene_id not in gene_sets:
                    gene_sets[gene_id] = set([gene_id])
            
            active_genes.add(gene_id)
        else:  # 'end'
            active_genes.discard(gene_id)

# Build overlap groups from gene_sets (extract unique sets)
print("Building overlap groups...", file=sys.stderr)
seen_sets = []
overlap_groups = []

for gene_id, gene_set in gene_sets.items():
    # Check if we've already recorded this set
    frozen_set = frozenset(gene_set)
    if frozen_set not in seen_sets:
        seen_sets.append(frozen_set)
        if len(gene_set) > 1:  # Only record actual overlaps
            overlap_groups.append(gene_set)

print(f"Found {len(overlap_groups)} overlap groups", file=sys.stderr)

# Apply filtering logic
print("\nApplying filtering logic...", file=sys.stderr)
genes_to_remove = set()

for group in overlap_groups:
    # Get gene info for all genes in group
    gene_info = []
    for gene_id in group:
        biotype = gene_biotypes.get(gene_id, "unknown")
        exon_count = exon_counts.get(gene_id, 0)
        gene_info.append((gene_id, biotype, exon_count))
    
    # Separate protein-coding from others
    protein_coding = [(g, b, e) for g, b, e in gene_info if b == "protein_coding"]
    non_coding = [(g, b, e) for g, b, e in gene_info if b != "protein_coding"]
    
    if protein_coding:
        # If there are protein-coding genes, remove all non-coding genes
        for gene_id, biotype, exon_count in non_coding:
            genes_to_remove.add(gene_id)
        
        # Among protein-coding genes, keep the one with most exons
        if len(protein_coding) > 1:
            protein_coding.sort(key=lambda x: x[2], reverse=True)  # Sort by exon count descending
            for gene_id, biotype, exon_count in protein_coding[1:]:
                genes_to_remove.add(gene_id)
    else:
        # All non-coding: keep the one with most exons
        if len(non_coding) > 1:
            non_coding.sort(key=lambda x: x[2], reverse=True)
            for gene_id, biotype, exon_count in non_coding[1:]:
                genes_to_remove.add(gene_id)

print(f"\nGenes to remove: {len(genes_to_remove)}", file=sys.stderr)

# Count what we're removing by biotype
remove_biotypes = collections.Counter([gene_biotypes.get(g, "unknown") for g in genes_to_remove])
print("\nRemoving by biotype:", file=sys.stderr)
for biotype, count in remove_biotypes.most_common():
    print(f"  {biotype}: {count}", file=sys.stderr)

# Separate into overlapping vs non-overlapping genes
overlapping_genes_all = set()
for group in overlap_groups:
    overlapping_genes_all.update(group)

genes_to_keep = set(gene_coords.keys()) - genes_to_remove
overlapping_kept = genes_to_keep & overlapping_genes_all

print(f"\nTotal overlapping genes: {len(overlapping_genes_all)}", file=sys.stderr)
print(f"Overlapping genes to keep: {len(overlapping_kept)}", file=sys.stderr)
print(f"Overlapping genes to remove: {len(genes_to_remove)}", file=sys.stderr)

# Write output files in the same directory as input GFF file
input_dir = os.path.dirname(input_gff)
if input_dir:
    keep_file = os.path.join(input_dir, f'{output_prefix}_to_keep.txt')
    remove_file = os.path.join(input_dir, f'{output_prefix}_to_remove.txt')
else:
    # Input file is in current directory
    keep_file = f'{output_prefix}_to_keep.txt'
    remove_file = f'{output_prefix}_to_remove.txt'

print(f"\nWriting output files...", file=sys.stderr)

# Genes to keep (only those that were in overlap groups)
with open(keep_file, 'w') as out:
    out.write("gene_id\tbiotype\texon_count\n")
    for gene_id in sorted(overlapping_kept):
        biotype = gene_biotypes.get(gene_id, "unknown")
        exon_count = exon_counts.get(gene_id, 0)
        out.write(f"{gene_id}\t{biotype}\t{exon_count}\n")

print(f"Written to {keep_file}: {len(overlapping_kept)} genes", file=sys.stderr)

# Genes to remove
with open(remove_file, 'w') as out:
    out.write("gene_id\tbiotype\texon_count\n")
    for gene_id in sorted(genes_to_remove):
        biotype = gene_biotypes.get(gene_id, "unknown")
        exon_count = exon_counts.get(gene_id, 0)
        out.write(f"{gene_id}\t{biotype}\t{exon_count}\n")

print(f"Written to {remove_file}: {len(genes_to_remove)} genes", file=sys.stderr)

# Count kept genes by biotype
keep_biotypes = collections.Counter([gene_biotypes.get(g, "unknown") for g in overlapping_kept])

# Generate statistics
print("\n=== SUMMARY ===", file=sys.stderr)
print(f"Total genes analyzed: {len(gene_coords)}", file=sys.stderr)
print(f"Genes with overlaps: {len(overlapping_genes_all)}", file=sys.stderr)
print(f"Overlap groups: {len(overlap_groups)}", file=sys.stderr)
print(f"Overlapping genes to keep: {len(overlapping_kept)}", file=sys.stderr)
print(f"Overlapping genes to remove: {len(genes_to_remove)}", file=sys.stderr)

print("\nKept overlapping genes by biotype:", file=sys.stderr)
for biotype, count in keep_biotypes.most_common(10):
    print(f"  {biotype}: {count}", file=sys.stderr)

print("\nRemoved genes by biotype:", file=sys.stderr)
for biotype, count in remove_biotypes.most_common(10):
    print(f"  {biotype}: {count}", file=sys.stderr)
