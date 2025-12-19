
import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt

# Define column names for BLAST tabular output (outfmt 6)
columns = ['qseqid', 'sseqid', '%identity', 'alignment_length', 'mismatches', 'gap_opens',
           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# Path to your BLAST result files
file_pattern = "./*.tblastn"  # Adjust this pattern to match your files

# Define thresholds
identity_threshold = 40.0  # Example: 90% identity
coverage_threshold = 0.80   # Example: 80% coverage (0.8)

# Initialize a dictionary to hold gene-genome relationships
gene_genome_dict = {}

# Get a list of all BLAST result files
blast_files = glob.glob(file_pattern)

# Read each BLAST result file
for file in blast_files:
    # Use the full file name (with extension) as the genome/database name
    genome_name = os.path.basename(file)  # Keeps the full file name with extension
    
    # Read the BLAST results
    blast_results = pd.read_csv(file, sep='\t', names=columns)
    
    # Loop through the BLAST results and filter by identity and coverage
    for _, row in blast_results.iterrows():
        gene = row['qseqid']  # Gene name (query sequence ID)

        # Calculate coverage
        query_length = row['qend'] - row['qstart'] + 1  # Length of the query sequence
        coverage = row['alignment_length'] / query_length  # Calculate coverage

        # Apply thresholds
        if row['%identity'] >= identity_threshold and coverage >= coverage_threshold:
            # Initialize dictionary for the gene if it doesn't exist
            if gene not in gene_genome_dict:
                gene_genome_dict[gene] = {}

            # Mark the genome (using filename with extension) as having a hit for the gene
            gene_genome_dict[gene][genome_name] = 1

# Get the full set of genes (rows) and genomes (columns)
all_genes = sorted(gene_genome_dict.keys())
all_genomes = sorted(set(genome for gene in gene_genome_dict for genome in gene_genome_dict[gene]))

# Create an empty dataframe with genomes as rows and genes as columns
presence_absence_matrix = pd.DataFrame(0, index=all_genomes, columns=all_genes)

# Populate the matrix with 1s where the gene is present in the genome
for gene, genomes in gene_genome_dict.items():
    for genome in genomes:
        presence_absence_matrix.at[genome, gene] = 1  # Transposed assignment

# Print the transposed presence/absence matrix
print(presence_absence_matrix)

# Optionally, save the matrix to a CSV file
presence_absence_matrix.to_csv('gene_presence_absence_transposed.csv')

# Create a heatmap of the transposed presence/absence matrix
plt.figure(figsize=(10, 8))  # Adjust the figure size as needed
sns.heatmap(presence_absence_matrix, cmap='viridis', cbar=True, linewidths=0.5)
plt.title('Gene Presence/Absence Heatmap (Transposed)')
plt.xlabel('Genes')
plt.ylabel('Genomes (TBLASTN Output Files)')
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
plt.tight_layout()  # Adjust layout to prevent clipping of tick-labels
plt.show()  # Display the heatmap

# Optionally, save the heatmap as an image file
plt.savefig('gene_presence_absence_heatmap_transposed.png', bbox_inches='tight')
