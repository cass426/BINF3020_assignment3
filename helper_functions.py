import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    Implement the `jukes_cantor` function in `helper_functions.py` to calculate the genetic distance 
    using Hamming distance and Jukes-Cantor correction. 
    Note that Jukes-Cantor correction ignores insertions and deletions, 
    so when calculating the length of the sequence, you should ignore positions with indels. 
    As a hint, JC corrected and non-corrected genetic distance (ignoring indels) 
    for SARS-CoV-2 sequences should be fairly similar and have the same unit.
    Example:
    ```
    AGGT-GA
    AGC-AGA
    ```
    This alignment has non-corrected genetic distance (ignoring indels) equal to 1, 
    while JC corrected genetic distance is 1.163.
    """
    valid_bases = set("ACGT")
    hamming_distance = 0
    valid_length = 0

    for ref_base, dist_base in zip(reference_sequence.upper(), distant_sequence.upper()):
        if ref_base not in valid_bases or dist_base not in valid_bases:
            continue
        valid_length += 1
        if ref_base != dist_base:
            hamming_distance += 1

    if valid_length == 0:
        return 0.0

    p = hamming_distance / valid_length

    if p >= 0.75:
        return float('inf')

    d_per_site = -0.75 * math.log(1 - (4 / 3) * p)
    return d_per_site * valid_length

print(jukes_cantor("AGGT-GA", "AGC-AGA"))

def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
