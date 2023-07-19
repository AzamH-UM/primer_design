# Primer Design

Utilities for creating primers for site directed mutagenesis

[Colab](https://colab.research.google.com/drive/1yuNZ3a_KaWDq-QeDM-aZHxYUHNFDqlD4?usp=sharing)

## Site Directed Mutagenesis Primer Design
primer_design.sdm.design_primer:

```
Design a forward primer for SDM.
    
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
    primer_df (pd.DataFrame): Dataframe containing all possible primers.
        
```

## Additional Utilities:

### Codon Table:
primer_design.codon_table
 - primer_design.codon_table['CODON'] -> 'AMINO ACID'
 
primer_design.codon_freq
 - primer_design.codon_freq['CODON'] -> FREQ
 
primer_design.highest_freq_codon(target_resn):
 - Get the highest frequency codon for a given amino acid
 - Args: target_resn (str): target amino acid
 - Returns: codon (str)
 - Example: highest_freq_codon('A') -> 'GCG'
 

### GC Content
primer_design.gc_content
 - primer_design.gc_content('GCAA') -> 0.5

### TM Calculator
primer_design.calculate_tm_nearest_neighbor(primer_sequence, primer_conc (in nM), polymerase)
SUPPORTED_POLYMERASES:   
 - Phusion
 - Phusion Hot Start Flex
 - Q5
 - Q5 Hot Start
 - Q5 Blood Direct
 - Q5U Hot Start
 - OneTaq
 - OneTaq Hot Start
 - Hot Start Taq
 - Taq
 - LongAmp Taq
 - LongAmp Hot Start Taq
 - Hemo KlenTaq
 - EpiMark Hot Start                        
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         



