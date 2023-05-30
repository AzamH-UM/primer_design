"""Utilities for sequence analysis"""


# Function to calculate GC content
def gc_content(primer):
    gc_count = primer.count("G") + primer.count("C")
    return gc_count / len(primer)
