"""Compute melting temperature of a primer sequence."""
import math
from primer_design.seq_utils import gc_content

_SUPPORTED_POLYMERASES = [
    "Phusion",
    "Phusion Hot Start Flex",
    "Q5",
    "Q5 Hot Start",
    "Q5 Blood Direct",
    "Q5U Hot Start",
    "OneTaq",
    "OneTaq Hot Start",
    "Hot Start Taq",
    "Taq",
    "LongAmp Taq",
    "LongAmp Hot Start Taq",
    "Hemo KlenTaq",
    "EpiMark Hot Start",
]


# Tm calculation function using nearest-neighbor method
def calculate_tm_nearest_neighbor(primer_sequence, primer_conc, polymerase):
    """Calculate melting temperature of a primer sequence using the nearest-neighbor method.

    Args:
        primer_sequence (str): Primer sequence.
        primer_conc (float): Primer concentration in M.
        polymerase (str): Polymerase to use.

    Returns:
        tm (float): Melting temperature in degrees Celsius."""

    # Check that polymerase is supported
    assert (
        polymerase in _SUPPORTED_POLYMERASES
    ), f"ERROR: Polymerase {polymerase} is not supported!"

    primer_sequence = primer_sequence.upper()

    enthalpy_dict = {
        "AA": -7.9,
        "TT": -7.9,
        "AT": -7.2,
        "TA": -7.2,
        "CA": -8.5,
        "TG": -8.5,
        "GT": -8.4,
        "AC": -8.4,
        "CT": -7.8,
        "AG": -7.8,
        "GA": -8.2,
        "TC": -8.2,
        "CG": -10.6,
        "GC": -9.8,
        "GG": -8.0,
        "CC": -8.0,
    }
    entropy_dict = {
        "AA": -22.2,
        "TT": -22.2,
        "AT": -20.4,
        "TA": -21.3,
        "CA": -22.7,
        "TG": -22.7,
        "GT": -22.4,
        "AC": -22.4,
        "CT": -21.0,
        "AG": -21.0,
        "GA": -22.2,
        "TC": -22.2,
        "CG": -27.2,
        "GC": -24.4,
        "GG": -19.9,
        "CC": -19.9,
    }

    enthalpy = 0
    entropy = 0

    for i in range(len(primer_sequence) - 1):
        pair = primer_sequence[i : i + 2]
        enthalpy += enthalpy_dict[pair]
        entropy += entropy_dict[pair]

    if primer_sequence[0] == "A" or primer_sequence[0] == "T":
        enthalpy += 2.3
        entropy += 4.1
    else:
        enthalpy += 0.1
        entropy += -2.8

    if primer_sequence[-1] == "A" or primer_sequence[-1] == "T":
        enthalpy += 2.3
        entropy += 4.1
    else:
        enthalpy += 0.1
        entropy += -2.8

    Tm_base = (1000 * enthalpy) / (entropy + 1.987 * math.log(primer_conc))

    cation_concentrations = {
        "Phusion": 160 / 1e3,
        "Phusion Hot Start Flex": 160 / 1e3,
        "Q5": 200 / 1e3,
        "Q5 Hot Start": 200 / 1e3,
        "Q5 Blood Direct": 200 / 1e3,
        "Q5U Hot Start": 250 / 1e3,
        "OneTaq": 55 / 1e3,
        "OneTaq Hot Start": 55 / 1e3,
        "Hot Start Taq": 55 / 1e3,
        "Taq": 55 / 1e3,
        "LongAmp Taq": 130 / 1e3,
        "LongAmp Hot Start Taq": 130 / 1e3,
        "Hemo KlenTaq": 70 / 1e3,
        "EpiMark Hot Start": 45 / 1e3,
    }

    m = cation_concentrations[polymerase]

    if polymerase == "Phusion" or polymerase == "Phusion Hot Start Flex":
        melting_temp = Tm_base - 273.15 + 16.6 * math.log10(m)
    else:
        melting_temp = (
            1
            / (
                (1 / Tm_base)
                + (
                    (4.29 * gc_content(primer_sequence) - 3.95) * math.log(m)
                    + 0.94 * (math.log(m)) ** 2
                )
                * 1e-05
            )
            - 273.15
        )

    return melting_temp


def calculate_tm_multi_sdm_agilent(primer_sequence, original_seq):
    """Calculate the melting temperature for multi-SDM"""
    gc_percent = gc_content(primer_sequence) * 100
    mismatch_percent = 100 - (
        100
        * sum([1 for i, j in zip(primer_sequence, original_seq) if i == j])
        / len(primer_sequence)
    )
    multi_sdm_tm = (
        81.5 + (0.41 * gc_percent) - (675 / len(primer_sequence)) - mismatch_percent
    )
    return multi_sdm_tm
