from Bio import SeqIO

# File paths
input_fasta = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results/intermediate_/filtered_no_outgroup_pruned.min4.fasta"
sample_file = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics/individuals_WCN.txt"
output_fasta = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Harmonia_genomics_2/results/intermediate_/filtered_samples.fasta"

# Read the list of sample IDs to keep
with open(sample_file, "r") as f:
    samples_to_keep = set(line.strip() for line in f)

# Filter FASTA records
filtered_records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    if record.id in samples_to_keep:
        filtered_records.append(record)

# Write the filtered records to a new FASTA file
SeqIO.write(filtered_records, output_fasta, "fasta")
print(f"Filtered FASTA saved to: {output_fasta}")
