#Load codon table
from primer_design.ecoli_codon_table import codon_table, codon_freq, highest_freq_codon

# Function to calculate GC content
from primer_design.seq_utils import gc_content

# Tm calculation function using nearest-neighbor method
from primer_design.tm_calculator import calculate_tm_nearest_neighbor

#Load site-directed mutagenesis package
from primer_design import sdm