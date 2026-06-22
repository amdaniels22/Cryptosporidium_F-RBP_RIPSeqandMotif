import pandas as pd
import numpy as np

# 1. Load the Master Table (entiched_genes_counts.txt)
# This defines the list of genes for the final output.
# Columns: gene_name, motif_counts
df_master = pd.read_csv('entiched_genes_counts.txt', sep='\t', header=None, names=['gene_name', 'motif_counts'])

# 2. Load Gene Lengths (all_mrna_cds_length.txt)
# Columns: gene_name, gene_length
df_len = pd.read_csv('all_mrna_cds_length.txt', sep='\t', header=None, names=['gene_name', 'gene_length'])

# 3. Load Input Gene List (input_gene_name.txt)
# Columns: gene_name
df_list = pd.read_csv('input_gene_name.txt', header=None, names=['gene_name'])

# 4. Load FIMO Data (fimo_cleaned.tsv)
# This file has a header.
df_fimo = pd.read_csv('fimo_cleaned.tsv', sep='\t')

# --- DATA PROCESSING ---

# Step A: Merge Lengths into Master Table
# We use a 'left' merge to keep only genes present in df_master
df_final = pd.merge(df_master, df_len, on='gene_name', how='left')

# Step B: Calculate Normalized Motif Counts (Column 4)
# Formula: motif_counts / gene_length per 1kb of gene
df_final['motif_counts_normalized_by_gene_length'] = (df_final['motif_counts'] / df_final['gene_length']) * 1000

# Step C: Determine Presence in List (Column 7)
# Check if gene_name is in the df_list dataframe
genes_in_list = set(df_list['gene_name'])
df_final['present_in_list'] = df_final['gene_name'].apply(lambda x: 'T' if x in genes_in_list else 'F')

# Step D: Process FIMO Data for Averages (Columns 5, 6, 8)
# We need to calculate statistics per gene from the FIMO file.

# First, merge gene length into FIMO data to calculate relative position
# We filter FIMO data to only keep genes present in our master list for efficiency
df_fimo_filtered = df_fimo[df_fimo['gene_name'].isin(df_master['gene_name'])].copy()
df_fimo_filtered = pd.merge(df_fimo_filtered, df_len, on='gene_name', how='left')

# Calculate Relative Position for every motif site
# Relative Position % = ((start + stop) / 2) / gene_length * 100
df_fimo_filtered['rel_position'] = ((df_fimo_filtered['start'] + df_fimo_filtered['stop']) / 2) / df_fimo_filtered['gene_length'] * 100

# Aggregate FIMO data by gene_name
# Calculate mean of p-value, q-value, and the derived relative position
fimo_stats = df_fimo_filtered.groupby('gene_name').agg(
    avg_p_value=('p-value', 'mean'),
    avg_q_value=('q-value', 'mean'),
    avg_motif_position=('rel_position', 'mean')
).reset_index()

# Step E: Merge FIMO Stats back into Master Table
df_final = pd.merge(df_final, fimo_stats, on='gene_name', how='left')

# --- FINAL OUTPUT FORMATTING ---

# Reorder columns to match user request
column_order = [
    'gene_name', 
    'motif_counts', 
    'gene_length', 
    'motif_counts_normalized_by_gene_length', 
    'avg_p_value', 
    'avg_motif_position', 
    'present_in_list', 
    'avg_q_value'
]

df_final = df_final[column_order]

# Save to file
output_filename = 'master_data_matrix.tsv'
df_final.to_csv(output_filename, sep='\t', index=False)

print(f"Processing complete. Output saved to {output_filename}")
print("\nPreview of the first 5 rows:")
print(df_final.head())