"""E coli Codon Table"""
import numpy as np

# Dictionary of codons and their corresponding amino acids
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

codon_freq = {
    'TTT': 1.9, 'TTC': 1.8, 'TTA': 1.0, 'TTG': 1.1,
    'TCT': 1.1, 'TCC': 1.0, 'TCA': 0.7, 'TCG': 0.8,
    'TAT': 1.6, 'TAC': 1.4, 'TAA': 0.2, 'TAG': .03,
    'TGT': 0.4, 'TGC': 0.6, 'TGA': 0.1, 'TGG': 1.4,

    'CTT': 1.0, 'CTC': 0.9, 'CTA': 0.3, 'CTG': 5.2,
    'CCT': 0.7, 'CCC': 0.4, 'CCA': 0.8, 'CCG': 2.4,
    'CAT': 1.2, 'CAC': 1.1, 'CAA': 1.3, 'CAG': 2.9,
    'CGT': 2.4, 'CGC': 2.2, 'CGA': 0.3, 'CGG': 0.5,

    'ATT': 2.7, 'ATC': 2.7, 'ATA': 0.4, 'ATG': 2.6,
    'ACT': 1.2, 'ACC': 2.4, 'ACA': 0.1, 'ACG': 1.3,
    'AAT': 1.6, 'AAC': 2.6, 'AAA': 3.8, 'AAG': 1.2,
    'AGT': 0.7, 'AGC': 1.5, 'AGA': 0.2, 'AGG': 0.2,

    'GTT': 2.0, 'GTC': 1.4, 'GTA': 1.2, 'GTG': 2.4,
    'GCT': 1.8, 'GCC': 2.3, 'GCA': 2.1, 'GCG': 3.2,
    'GAT': 3.3, 'GAC': 2.3, 'GAA': 4.4, 'GAG': 1.9,
    'GGT': 2.8, 'GGC': 3.0, 'GGA': 0.7, 'GGG': 0.9,
}

def highest_freq_codon(target_resn):
    """Get the highest frequency codon for a given amino acid
    
    Args:
        target_resn (str): target amino acid
        
    Returns:
        str: codon
    """
    target_codons = [codon for codon, resn in codon_table.items() if resn == target_resn]
    target_freqs = [codon_freq[codon] for codon in target_codons]
    return target_codons[np.argmax(target_freqs)]
    