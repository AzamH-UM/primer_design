"""Perform site-directed mutagenesis (SDM) on a DNA sequence."""
import primer_design 
import pandas as pd 
from Bio.Seq import Seq
from copy import deepcopy

def design_primer(dna_seq,
                          wt_resn,
                          mut_resi,
                          mut_resn,
                          max_pad = 20,
                          min_pad = 10,
                          primer_conc = 500,
                          GC_3prime = True,
                          GC_5prime = False,
                          TM_target = 70,
                          polymerase = 'Phusion',
                          verbose = False,
                          GC_max = 0.8,
                          TM_diff_max = None,
                          show_gap = True,
                          ):

    """Design a forward primer for SDM.
    
    Args:
        dna_seq (str): DNA sequence to mutate.
        wt_resn (str): Wildtype residue.
        mut_resi (int): Residue position to mutate.
        mut_resn (str): Mutant residue.
        max_pad (int): Maximum number of nucleotides to pad the codon with. Default is 20.
        min_pad (int): Minimum number of nucleotides to pad the codon with. Default is 10.
        primer_conc (float): Primer concentration in nM. Default is 500.
        GC_3prime (bool): Filter forward primers so that 3' is G or C . Default is True.
        GC_5prime (bool): Filter forward primers so that 5' is G or C . Default is False.
        TM_target (float): Target melting temperature. Default is 70.
        polymerase (str): Polymerase used for PCR. Default is 'Phusion'.
        verbose (bool): Print information about primer design. Default is False.
        GC_max (float): Maximum GC content. Default is 0.8.
        TM_diff_max (float): Optional: Maximum allowable difference in melting temperature between forward and reverse primers. Default is None.
        show_gap (bool): Add gap to primer to indicate where codon is. Default is True.
        
    Returns:
        primer_df (pd.DataFrame): Dataframe containing all possible primers."""
    
    primer_conc = primer_conc / 1e9 #Convert to M
    
    #Find dna codon in sequence
    codon_start, codon_end = (mut_resi - 1) * 3, (mut_resi - 1) * 3 + 3
    wt_codon = dna_seq[codon_start:codon_end]
    if not primer_design.codon_table[wt_codon] == wt_resn:
        print(f'{mut_resi=}')
        print(f'{dna_seq=}')
        print(f'{wt_codon=}')
        print(f'{wt_resn=}')
        print(f'{primer_design.codon_table[wt_codon]=}')
        raise ValueError('ERROR: Codon at residue position does not translate to wildtype residue!')

    #Get codon for mutation
    mut_codon = primer_design.highest_freq_codon(mut_resn)

    #Mutate codon
    original_seq = deepcopy(dna_seq)
    dna_seq = dna_seq[:codon_start] + mut_codon + dna_seq[codon_end:]
    
    #Create forward primer by slicing max-pad - codon - max-pad
    forward_primer = dna_seq[codon_start - max_pad: codon_end + max_pad]

    #Create dataframe to store all primers
    columns = ['Forward', 'Reverse', 'TM', 'Reverse TM', 'GC %',]

    #Filter columns:
    filter_columns = ['Ends G/C', 'Starts G/C',
                      'TM Difference',]
    primer_df = pd.DataFrame(columns = columns + filter_columns)
    
    #Loop through all possible primers and classify as candidates
    for five_prime_cut in range(max_pad - min_pad + 1):
        for three_prime_cut in range(max_pad - min_pad + 1):
            primer_name = f'{max_pad - five_prime_cut}-codon-{max_pad - three_prime_cut}'
            primer = forward_primer[five_prime_cut: len(forward_primer) - three_prime_cut]
            reverse_primer = str(Seq(primer).reverse_complement())
            tm = primer_design.calculate_tm_nearest_neighbor(primer, primer_conc, polymerase)
            tm_reverse = primer_design.calculate_tm_nearest_neighbor(reverse_primer, primer_conc, polymerase)
            gc = primer_design.gc_content(primer)

            #Insert gap to show where codon is
            if show_gap:
                primer_gapped = list(primer)
                primer_gapped.insert(max_pad - five_prime_cut, '-')
                primer_gapped.insert(max_pad - five_prime_cut + 4, '-')
                primer_gapped = ''.join(primer_gapped)

                reverse_gapped = list(reverse_primer)
                reverse_gapped.insert(max_pad - three_prime_cut, '-')
                reverse_gapped.insert(max_pad - three_prime_cut + 4, '-')
                reverse_gapped = ''.join(reverse_gapped)
            else:
                primer_gapped = primer
                reverse_gapped = reverse_primer

            #Add to dataframe
            primer_df.at[primer_name, 'Forward'] = primer_gapped
            primer_df.at[primer_name, 'Reverse'] = reverse_gapped
            primer_df.at[primer_name, 'TM'] = tm
            primer_df.at[primer_name, 'Reverse TM'] = tm_reverse
            primer_df.at[primer_name, 'GC %'] = gc

            #Check if 5' and 3' ends are G or C
            ends_gc = False
            starts_gc = False

            if primer[0] == 'G' or primer[0] == 'C':
                starts_gc = True
            if primer[-1] == 'G' or primer[-1] == 'C':
                ends_gc = True
            primer_df.at[primer_name, 'Ends G/C'] = ends_gc
            primer_df.at[primer_name, 'Starts G/C'] = starts_gc

            #Calculate TM difference
            tm_diff = abs(tm - tm_reverse)
            primer_df.at[primer_name, 'TM Difference'] = tm_diff
            
    
    #Filter primers by TM and GC content
    if GC_5prime:
        primer_df = primer_df[primer_df['Starts G/C'] == True]
    if GC_3prime:
        primer_df = primer_df[primer_df['Ends G/C'] == True]

    #Filter by TM difference
    if TM_diff_max:
        primer_df = primer_df[primer_df['TM Difference'] <= TM_diff_max]


    #Filter by GC content
    primer_df = primer_df[primer_df['GC %'] <= GC_max]

    #Sort by closest TM to target
    primer_df = primer_df.iloc[(primer_df['TM']-TM_target).abs().argsort()]

    if verbose:
        print(f'{original_seq=}')
        print(f'{wt_resn=}, {mut_resi=}, {mut_resn=}')
        print(f'{wt_codon=}')
        print(f'{mut_codon=}')

    return primer_df[columns]