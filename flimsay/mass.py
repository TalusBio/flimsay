from .constants import PROTON, STD_AA_MASS, WATER


def mz_to_mass(mz, charge):
    """Calculate the mass of a peptide given its m/z and charge.

    Examples
    --------
    >>> mz_to_mass(500.000, 2)
    997.98434993586
    """
    return (mz * charge) - charge * PROTON


def mass_to_mz(mass, charge):
    """Calculate the m/z of a peptide given its mass and charge.

    This assumes the charge adduct is a proton.

    Examples
    --------
    >>> mass_to_mz(1000.000, 2)
    501.00782503207
    """
    return (mass + charge * PROTON) / charge


# ACLLPEPTK
# 1027.5372 = Mass for carbamidomethylated version of the peptide
# From https://web.expasy.org/peptide_mass/


def calc_mass(sequence) -> float:
    """Calculate the mass of a peptide sequence.

    Parameters
    ----------
    sequence : str
        Peptide

    Returns
    -------
    float
        Mass of the peptide, no protons and assuming
        carbamidomethylation of cysteine residues.

    Examples
    --------
    >>> calc_mass("ACLLPEPTK")
    1027.53721747546
    """
    mass = WATER
    for aa in sequence:
        try:
            mass += STD_AA_MASS[aa]
        except KeyError:
            raise KeyError(f"Unknown amino acid: {aa} in peptide: {sequence}")
    return mass
