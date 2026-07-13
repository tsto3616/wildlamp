from Bio import SeqIO
import math

# Parse the FASTA file and get the lengths of all sequences
lengths = [len(record.seq) for record in SeqIO.parse("HSP70_all.fasta", "fasta")]

# Calculate the required metrics
total_base_pairs = sum(lengths)
num_seqs = len(lengths)

mean_length = total_base_pairs / num_seqs if num_seqs > 0 else 0
min_length = min(lengths) if lengths else 0
max_length = max(lengths) if lengths else 0

# Standard deviation
if num_seqs > 1:
    variance = sum((L - mean_length)**2 for L in lengths) / (num_seqs - 1)
    std_dev = math.sqrt(variance)
else:
    std_dev = 0

# Print the results
print(f"Total Base Pairs: {total_base_pairs}")
print(f"Mean (Average) Length: {mean_length:.2f}")
print(f"Standard Deviation: {std_dev:.2f}")
print(f"Range: {min_length} to {max_length}")
