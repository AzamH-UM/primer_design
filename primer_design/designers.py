"""Perform site-directed mutagenesis (SDM) on a DNA sequence."""
import primer_design
import pandas as pd
from Bio.Seq import Seq
from copy import deepcopy


def _resi_to_codon(dna_seq, resi):
    codon_start, codon_end = (resi - 1) * 3, (resi - 1) * 3 + 3
    codon = dna_seq[codon_start:codon_end]
    return codon


def sdm(
    dna_seq,
    wt_resn,
    mut_resi,
    mut_resn,
    max_pad=20,
    min_pad=10,
    primer_conc=500,
    TM_target=70,
    GC_min=0.0,
    GC_max=0.8,
    GC_3prime=True,
    GC_5prime=False,
    polymerase="Phusion",
    TM_diff_max=None,
    show_gap=True,
    verbose=False,
):
    """Design forward and reverse primer for Site-Directed Mutagenesis (SDM).

    Args:
        dna_seq (str): DNA sequence to mutate.
        wt_resn (str): Wildtype residue.
        mut_resi (int): Residue position to mutate.
        mut_resn (str): Mutant residue.
        max_pad (int): Maximum number of nucleotides to pad the codon with. Default is 20.
        min_pad (int): Minimum number of nucleotides to pad the codon with. Default is 10.
        primer_conc (float): Primer concentration in nM. Default is 500.
        TM_target (float): Target melting temperature. Default is 70.
        GC_min (float): Minimum GC content. Default is 0.0.
        GC_max (float): Maximum GC content. Default is 0.8.
        GC_3prime (bool): Filter forward primers so that 3' is G or C . Default is True.
        GC_5prime (bool): Filter forward primers so that 5' is G or C . Default is False.
        polymerase (str): Polymerase used for PCR. Default is 'Phusion'.
        TM_diff_max (float): Optional: Maximum allowable difference in melting temperature between forward and reverse primers. Default is None.
        show_gap (bool): Add gap to primer to indicate where codon is. Default is True.
        verbose (bool): Print information about primer design. Default is False.

    Returns:
        primer_df (pd.DataFrame): Dataframe containing all possible primers.
    """
    primer_conc = primer_conc / 1e9  # Convert to M

    # Find dna codon in sequence
    codon_start, codon_end = (mut_resi - 1) * 3, (mut_resi - 1) * 3 + 3
    wt_codon = dna_seq[codon_start:codon_end]
    if not primer_design.codon_table[wt_codon] == wt_resn:
        print(f"{mut_resi=}")
        print(f"{dna_seq=}")
        print(f"{wt_codon=}")
        print(f"{wt_resn=}")
        print(f"{primer_design.codon_table[wt_codon]=}")
        raise ValueError(
            "ERROR: Codon at residue position does not translate to wildtype residue!"
        )

    # Get codon for mutation
    mut_codon = primer_design.highest_freq_codon(mut_resn)

    # Mutate codon
    original_seq = deepcopy(dna_seq)
    dna_seq = dna_seq[:codon_start] + mut_codon + dna_seq[codon_end:]

    # Create forward primer by slicing max-pad - codon - max-pad
    forward_primer = dna_seq[codon_start - max_pad : codon_end + max_pad]

    # Create dataframe to store all primers
    columns = [
        "Forward",
        "Reverse",
        "Forward TM",
        "Reverse TM",
        "GC",
    ]

    # Filter columns:
    filter_columns = [
        "Ends G/C",
        "Starts G/C",
        "TM Difference",
    ]

    primer_df = pd.DataFrame(columns=columns + filter_columns)

    # Loop through all possible primers and classify as candidates
    for five_prime_cut in range(max_pad - min_pad + 1):
        for three_prime_cut in range(max_pad - min_pad + 1):
            # Name of primer is padding-codon-padding
            primer_name = (
                f"{max_pad - five_prime_cut}-codon-{max_pad - three_prime_cut}"
            )

            # Slice primer
            primer = forward_primer[
                five_prime_cut : len(forward_primer) - three_prime_cut
            ]
            reverse_primer = str(Seq(primer).reverse_complement())

            # Compute TM and GC content
            tm = primer_design.calculate_tm_nearest_neighbor(
                primer, primer_conc, polymerase
            )
            tm_reverse = primer_design.calculate_tm_nearest_neighbor(
                reverse_primer, primer_conc, polymerase
            )
            gc = primer_design.gc_content(primer)

            # Insert gap to show where codon is
            if show_gap:
                primer_gapped = list(primer)
                primer_gapped.insert(max_pad - five_prime_cut, "-")
                primer_gapped.insert(max_pad - five_prime_cut + 4, "-")
                primer_gapped = "".join(primer_gapped)

                reverse_gapped = list(reverse_primer)
                reverse_gapped.insert(max_pad - three_prime_cut, "-")
                reverse_gapped.insert(max_pad - three_prime_cut + 4, "-")
                reverse_gapped = "".join(reverse_gapped)
            else:
                primer_gapped = primer
                reverse_gapped = reverse_primer

            # Add to dataframe
            primer_df.at[primer_name, "Forward"] = primer_gapped
            primer_df.at[primer_name, "Reverse"] = reverse_gapped
            primer_df.at[primer_name, "Forward TM"] = tm
            primer_df.at[primer_name, "Reverse TM"] = tm_reverse
            primer_df.at[primer_name, "GC"] = gc

            # Check if 5' and 3' ends are G or C
            ends_gc = False
            starts_gc = False

            if primer[0] == "G" or primer[0] == "C":
                starts_gc = True
            if primer[-1] == "G" or primer[-1] == "C":
                ends_gc = True
            primer_df.at[primer_name, "Ends G/C"] = ends_gc
            primer_df.at[primer_name, "Starts G/C"] = starts_gc

            # Calculate TM difference
            tm_diff = abs(tm - tm_reverse)
            primer_df.at[primer_name, "TM Difference"] = tm_diff

    # Filter primers by TM and GC content
    if GC_5prime:
        primer_df = primer_df[primer_df["Starts G/C"] == True]
    if GC_3prime:
        primer_df = primer_df[primer_df["Ends G/C"] == True]

    # Filter by TM difference
    if TM_diff_max:
        primer_df = primer_df[primer_df["TM Difference"] <= TM_diff_max]

    # Filter by GC content
    primer_df = primer_df[primer_df["GC"] <= GC_max]
    primer_df = primer_df[primer_df["GC"] >= GC_min]

    # Sort by closest TM to target
    primer_df = primer_df.iloc[(primer_df["Forward TM"] - TM_target).abs().argsort()]

    if verbose:
        print(f"{original_seq=}")
        print(f"{wt_resn=}, {mut_resi=}, {mut_resn=}")
        print(f"{wt_codon=}")
        print(f"{mut_codon=}")

    return primer_df[columns]


def sdm_multi(
    dna_seq,
    wt_resn,
    mut_resi,
    mut_resn,
    max_pad=21,
    min_pad=11,
    TM_target=80,
    GC_min=0.4,
    GC_max=1.0,
    GC_3prime=True,
    GC_5prime=False,
    show_gap=True,
    verbose=False,
):
    """Design forward and reverse primer for multi-site SDM.

    Args:
        dna_seq (str): DNA sequence to mutate.
        wt_resn (str): Wildtype residue.
        mut_resi (int): Residue position to mutate.
        mut_resn (str): Mutant residue.
        max_pad (int): Maximum number of nucleotides to pad the codon with.
        min_pad (int): Minimum number of nucleotides to pad the codon with.
        TM_target (float): Target melting temperature.
        GC_min (float): Minimum GC content.
        GC_max (float): Maximum GC content.
        GC_3prime (bool): Whether to filter primers by 3' GC content.
        GC_5prime (bool): Whether to filter primers by 5' GC content.
        show_gap (bool): Whether to show gap in primer.
        verbose (bool): Whether to print primer design information.
    """

    # Get codon
    codon_start, codon_end = (mut_resi - 1) * 3, (mut_resi - 1) * 3 + 3
    wt_codon = dna_seq[codon_start:codon_end]
    mut_codon = primer_design.highest_freq_codon(mut_resn)

    # Get primer df from sdm
    primer_df = sdm(
        dna_seq=dna_seq,
        wt_resn=wt_resn,
        mut_resi=mut_resi,
        mut_resn=mut_resn,
        max_pad=max_pad,
        min_pad=min_pad,
        TM_target=TM_target,
        GC_min=GC_min,
        GC_max=1.0,
        GC_3prime=GC_3prime,
        GC_5prime=GC_5prime,
        show_gap=show_gap,
        verbose=verbose,
    )

    # Recompute TM

    for index, row in primer_df.iterrows():
        five_prime_len = int(index.split("-")[0])
        three_prime_len = int(index.split("-")[-1])
        if show_gap:
            five_prime_len += 1
            three_prime_len += 1
        forward_primer = row["Forward"]
        reverse_primer = row["Reverse"]
        original_forward = (
            forward_primer[0:five_prime_len]
            + wt_codon
            + forward_primer[-three_prime_len:]
        )
        original_reverse = str(Seq(original_forward).reverse_complement())
        tm = primer_design.calculate_tm_multi_sdm_agilent(
            forward_primer, original_forward
        )
        tm_reverse = primer_design.calculate_tm_multi_sdm_agilent(
            reverse_primer, original_reverse
        )
        primer_df.at[index, "Forward TM"] = tm
        primer_df.at[index, "Reverse TM"] = tm_reverse

    # Sort by closest TM to target
    primer_df = primer_df.iloc[(primer_df["Forward TM"] - TM_target).abs().argsort()]

    return primer_df


def ssm_single(
    dna_seq,
    wt_resn,
    mut_resi,
    max_pad=20,
    min_pad=10,
    primer_conc=500,
    TM_target=70,
    GC_min=0.0,
    GC_max=0.8,
    GC_3prime=True,
    GC_5prime=False,
    polymerase="Phusion",
    TM_diff_max=None,
    show_gap=True,
    verbose=False,
):
    """Generate single-site saturation mutagenesis single primer.

    Args:
        dna_seq (str): DNA sequence to mutate.
        wt_resn (str): Wildtype residue.
        mut_resi (int): Residue position to mutate.
        max_pad (int): Maximum number of nucleotides to pad the codon with.
        min_pad (int): Minimum number of nucleotides to pad the codon with.
        primer_conc (int): Primer concentration in nM.
        TM_target (float): Target melting temperature.
        GC_min (float): Minimum GC content.
        GC_max (float): Maximum GC content.
        GC_3prime (bool): Whether to filter primers by 3' GC content.
        GC_5prime (bool): Whether to filter primers by 5' GC content.
        polymerase (str): Polymerase to use for TM calculation.
        TM_diff_max (float): Maximum difference in TM between forward and reverse primers.
        show_gap (bool): Whether to show gap in primer.
        verbose (bool): Whether to print primer design information.

    """

    # Generate primers with sdm
    primer_df = sdm(
        dna_seq=dna_seq,
        wt_resn=wt_resn,
        mut_resi=mut_resi,
        mut_resn=wt_resn,
        max_pad=max_pad,
        min_pad=min_pad,
        TM_target=TM_target,
        GC_min=GC_min,
        GC_max=GC_max,
        GC_3prime=GC_3prime,
        GC_5prime=GC_5prime,
        show_gap=show_gap,
        verbose=verbose,
    )

    # Swap out wt codon for NNS codon
    mut_codon = "NNS"

    for index, row in primer_df.iterrows():
        five_prime_len = int(index.split("-")[0])
        three_prime_len = int(index.split("-")[-1])
        if show_gap:
            five_prime_len += 1
            three_prime_len += 1
        forward_primer = row["Forward"]
        reverse_primer = row["Reverse"]
        new_forward = (
            forward_primer[0:five_prime_len]
            + mut_codon
            + forward_primer[-three_prime_len:]
        )
        primer_df.at[index, "Forward"] = new_forward

    # Drop reverse primers
    primer_df = primer_df.drop(columns=["Reverse", "Reverse TM"])

    return primer_df


def ssm_double(
    dna_seq,
    wt_resn,
    mut_resi,
    max_len=30,
    min_len=10,
    TM_target=64,
    TM_range=2,
    GC_min=0.0,
    GC_max=1.0,
    GC_end=True,
    polymerase="Phusion",
    primer_conc=500,
):
    """Design forward and reverse primer for site-saturation mutagenesis

    Args:
        dna_seq (str): DNA sequence to mutate.
        wt_resn (str): Wildtype residue.
        mut_resi (int): Residue position to mutate.
        max_len (int): Maximum length of primer.
        min_len (int): Minimum length of primer.
        TM_target (float): Target melting temperature.
        GC_min (float): Minimum GC content.
        GC_max (float): Maximum GC content.
        GC_end (bool): Whether to force G/C at end
        polymerase (str): Polymerase to use for TM calculation.
        primer_conc (int): Primer concentration in nM.
        verbose (bool): Whether to print primer design information.

    """

    primer_conc = primer_conc / 1e9  # Convert to M

    # Create dataframe to store primers
    columns = [
        "Sequence",
        "Length",
        "GC Content",
        "TM",
    ]
    filter_columns = [
        "GC end",
        "Length",
        "GC within range",
        "TM within range",
    ]
    for_primer_df = pd.DataFrame(columns=columns + filter_columns)
    rev_primer_df = pd.DataFrame(columns=columns + filter_columns)

    # Get wildtype codon
    wt_codon = dna_seq[3 * (mut_resi - 1) : 3 * mut_resi]

    # Get degenerate codon
    degenerate_codon = "NNK"

    # Create forward primer
    forward_start = 3 * mut_resi
    forward_primer = dna_seq[
        forward_start : min(forward_start + (max_len - 3), len(dna_seq))
    ]

    # Create reverse primer
    reverse_primer = dna_seq[
        max(forward_start - 3 - (max_len - 3), 0) : forward_start - 3
    ]
    reverse_primer = "".join(str(Seq(reverse_primer).reverse_complement()))

    # Loop over possible primers
    for i_slice in range((max_len - min_len - 3)):
        forward_primer_candidate = forward_primer[0 : len(forward_primer) - i_slice]
        reverse_primer_candidate = reverse_primer[0 : len(reverse_primer) - i_slice]

        # Compute TM
        tm = primer_design.calculate_tm_nearest_neighbor(
            forward_primer_candidate, primer_conc, polymerase
        )
        tm_rev = primer_design.calculate_tm_nearest_neighbor(
            reverse_primer_candidate, primer_conc, polymerase
        )

        gc = primer_design.gc_content(forward_primer_candidate)
        gc_rev = primer_design.gc_content(reverse_primer_candidate)

        # Check if TM is within range
        tm_in_range = False
        rev_tm_in_range = False
        if TM_target - TM_range <= tm <= TM_target + TM_range:
            tm_in_range = True
        if TM_target - TM_range <= tm_rev <= TM_target + TM_range:
            rev_tm_in_range = True

        # Check if GC is within range
        gc_in_range = False
        rev_gc_in_range = False
        if GC_min <= gc <= GC_max:
            gc_in_range = True
        if GC_min <= gc_rev <= GC_max:
            rev_gc_in_range = True

        # Check that forward and reverse primers end in G/C
        gc_end = False
        gc_rev_end = False
        if forward_primer_candidate[-1] == "G" or forward_primer_candidate[-1] == "C":
            gc_end = True
        if reverse_primer_candidate[-1] == "G" or reverse_primer_candidate[-1] == "C":
            gc_rev_end = True

        # Add degenerate codon
        forward_primer_candidate = degenerate_codon + forward_primer_candidate

        # Add to dataframe
        for_primer_df.at[i_slice, "Sequence"] = forward_primer_candidate
        for_primer_df.at[i_slice, "Length"] = len(forward_primer_candidate)
        for_primer_df.at[i_slice, "GC Content"] = gc
        for_primer_df.at[i_slice, "TM"] = tm
        for_primer_df.at[i_slice, "GC end"] = gc_end
        for_primer_df.at[i_slice, "GC within range"] = gc_in_range
        for_primer_df.at[i_slice, "TM within range"] = tm_in_range

        rev_primer_df.at[i_slice, "Sequence"] = reverse_primer_candidate
        rev_primer_df.at[i_slice, "Length"] = len(reverse_primer_candidate)
        rev_primer_df.at[i_slice, "GC Content"] = gc_rev
        rev_primer_df.at[i_slice, "TM"] = tm_rev
        rev_primer_df.at[i_slice, "GC end"] = gc_rev_end
        rev_primer_df.at[i_slice, "GC within range"] = rev_gc_in_range
        rev_primer_df.at[i_slice, "TM within range"] = rev_tm_in_range

    # Filter primers
    for_primer_df = for_primer_df[
        (for_primer_df["GC end"] == GC_end)
        & (for_primer_df["GC within range"] == True)
        & (for_primer_df["TM within range"] == True)
    ]

    rev_primer_df = rev_primer_df[
        (rev_primer_df["GC end"] == GC_end)
        & (rev_primer_df["GC within range"] == True)
        & (rev_primer_df["TM within range"] == True)
    ]

    # Sort by closest to target TM
    for_primer_df = for_primer_df.sort_values(
        by=["TM"], key=lambda x: abs(x - TM_target), ascending=True
    )
    rev_primer_df = rev_primer_df.sort_values(
        by=["TM"], key=lambda x: abs(x - TM_target), ascending=True
    )

    # Remove filter columns
    for_primer_df = for_primer_df[columns]
    rev_primer_df = rev_primer_df[columns]

    return for_primer_df, rev_primer_df


# Design primers for a range of residues
def ssm_multi(
    dna_seq,
    residues,
    ssm_kwargs,
):
    """Design primers for a range of residues

    Args:
        dna_seq (str): DNA sequence
        residues (list): List of residues to mutate
        ssm_kwargs (dict): Keyword arguments for ssm_double

    Returns:
        primers (pd.DataFrame): List of primers

    """

    primers = []

    # Translate DNA sequence
    aa_seq = Seq(dna_seq).translate()

    # Iterate over residues
    for resi in residues:
        wt_resn = aa_seq[resi - 1]

        # Design primers
        for_primers, rev_primers = primer_design.ssm_double(
            dna_seq=dna_seq, wt_resn=wt_resn, mut_resi=resi, **ssm_kwargs
        )

        # Select first forward primer and reverse primer
        for_primer = None
        rev_primer = None
        for_TM = None
        rev_TM = None
        if not for_primers.empty:
            for_primer = for_primers.iloc[0]["Sequence"]
            for_TM = for_primers.iloc[0]["TM"]
        if not rev_primers.empty:
            rev_primer = rev_primers.iloc[0]["Sequence"]
            rev_TM = rev_primers.iloc[0]["TM"]

        # Convert primers to strings
        for_primer = str(for_primer)
        rev_primer = str(rev_primer)

        # Add to list
        primers.append((for_primer, rev_primer, for_TM, rev_TM))

    # Convert to dataframe
    primers = pd.DataFrame(primers, columns=["forward", "reverse", "for_tm", "rev_tm"], index=residues)
    primers.index.name = "residue"

    return primers


def insertion(
    dna_seq,
    ins_resi,
    ins_resn,
    max_len=30,
    min_len=10,
    TM_target=63,
    GC_min=0.0,
    GC_max=1.0,
    GC_3prime=True,
    GC_5prime=False,
    polymerase="Phusion",
    primer_conc=500,
    verbose=False,
    show_gap=True,
):
    """Design forward and reverse primer for insertion.

    Args:
        dna_seq (str): DNA sequence to insert.
        ins_resi (int): Residue position to insert residue before.
        ins_resn (str): Residue to insert.
        max_len (int): Maximum length of primer.
        min_len (int): Minimum length of primer.
        TM_target (float): Target melting temperature.
        GC_min (float): Minimum GC content.
        GC_max (float): Maximum GC content.
        GC_3prime (bool): Whether to filter primers by 3' GC content.
        GC_5prime (bool): Whether to filter primers by 5' GC content.
        polymerase (str): Polymerase to use for TM calculation.
        primer_conc (float): Primer concentration in nM.
        verbose (bool): Whether to print primer design information.
        show_gap (bool): Whether to show gap in primer.

    """

    primer_conc = primer_conc / 1e9  # Convert to M

    # Create dataframe to store all primers
    columns = [
        "Sequence",
        "TM",
        "GC",
    ]

    # Filter columns:
    filter_columns = [
        "Ends G/C",
        "Starts G/C",
    ]

    primer_df = pd.DataFrame(columns=columns + filter_columns)

    # Get insertion codon
    ins_codon = primer_design.highest_freq_codon(ins_resn)

    # Create forward primer
    forward_start = ins_resi * 3 - 3
    forward_primer = ins_codon + dna_seq[forward_start : forward_start + (max_len - 3)]

    # Create reverse primer
    reverse_primer = dna_seq[forward_start - (max_len - 3) : forward_start]

    # Explore all possible lengths
    for i_slice in range((max_len - min_len - 3)):
        forward_primer_candidate = forward_primer[0 : len(forward_primer) - i_slice]
        reverse_primer_candidate = reverse_primer[i_slice:]

        # Compute lengths
        forward_len = len(forward_primer_candidate)
        reverse_len = len(reverse_primer_candidate)

        # Create entry names
        forward_name = f"forward-{forward_len}"
        reverse_name = f"reverse-{reverse_len}"

        # Compute TM
        forward_tm = primer_design.calculate_tm_nearest_neighbor(
            forward_primer_candidate, primer_conc, polymerase
        )
        reverse_tm = primer_design.calculate_tm_nearest_neighbor(
            reverse_primer_candidate, primer_conc, polymerase
        )

        # Compute GC content
        forward_gc = primer_design.gc_content(forward_primer_candidate)
        reverse_gc = primer_design.gc_content(reverse_primer_candidate)

        # Determine if primers ends/starts with G/C
        forward_ends_gc = forward_primer_candidate[-1] in ["G", "C"]
        reverse_ends_gc = reverse_primer_candidate[-1] in ["G", "C"]
        forward_starts_gc = forward_primer_candidate[0] in ["G", "C"]
        reverse_starts_gc = reverse_primer_candidate[0] in ["G", "C"]

        # Add gap to forward primer
        if show_gap:
            forward_primer_candidate = (
                forward_primer_candidate[0:3] + "-" + forward_primer_candidate[3:]
            )

        # Add to dataframe
        primer_df.at[forward_name, "Sequence"] = forward_primer_candidate
        primer_df.at[forward_name, "TM"] = forward_tm
        primer_df.at[forward_name, "GC"] = forward_gc
        primer_df.at[reverse_name, "Sequence"] = reverse_primer_candidate
        primer_df.at[reverse_name, "TM"] = reverse_tm
        primer_df.at[reverse_name, "GC"] = reverse_gc
        primer_df.at[forward_name, "Ends G/C"] = forward_ends_gc
        primer_df.at[reverse_name, "Ends G/C"] = reverse_ends_gc
        primer_df.at[forward_name, "Starts G/C"] = forward_starts_gc
        primer_df.at[reverse_name, "Starts G/C"] = reverse_starts_gc

    # Filter by GC content
    primer_df = primer_df[primer_df["GC"] <= GC_max]
    primer_df = primer_df[primer_df["GC"] >= GC_min]

    # Filter primers by TM and GC content
    if GC_5prime:
        primer_df = primer_df[primer_df["Starts G/C"] == True]
    if GC_3prime:
        primer_df = primer_df[primer_df["Ends G/C"] == True]

    # Sort by closest TM to target
    primer_df = primer_df.iloc[(primer_df["TM"] - TM_target).abs().argsort()]

    if verbose:
        print(f"{dna_seq=}")
        print(f"{ins_resi=}, {ins_resn=}")
        print(f"{ins_codon=}")

    return primer_df[columns]
