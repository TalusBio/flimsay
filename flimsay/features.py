import re

import lightgbm as lgb
import numpy as np
from loguru import logger

from .constants import BULKY_AAS, NEGATIVE_AAS, POSITIVE_AAS
from .data import select_split
from .mass import calc_mass, mass_to_mz, mz_to_mass  # noqa

FEATURE_COLUMN_DESCRIPTIONS = {
    "PrecursorMz": "Measured precursor m/z",
    "Mass": "Measured precursor mass (Da)",
    "PrecursorCharge": "Measured precursor charge, from the isotope envelope",
    "PepLength": "Length of the peptide sequence in amino acids",
    "NumBulky": f"Number of bulky amino acids ({BULKY_AAS})",
    "NumPos": f"Number of positive amino acids ({POSITIVE_AAS})",
    "PosIndexL": f"Relative position of the first positive amino acid ({POSITIVE_AAS})",
    "PosIndexR": f"Relative position of the last positive amino acid ({POSITIVE_AAS})",
    "NumNeg": f"Number of negative amino acids ({NEGATIVE_AAS})",
    "NegIndexL": f"Relative position of the first negative amino acid ({NEGATIVE_AAS})",
    "NegIndexR": f"Relative position of the last negative amino acid ({NEGATIVE_AAS})",
}

FEATURE_COLUMNS = list(FEATURE_COLUMN_DESCRIPTIONS)

POSITIVE_AAS_PATTERN = f"[{POSITIVE_AAS}]"
NEGATIVE_AAS_PATTERN = f"[{NEGATIVE_AAS}]"
BULKY_AAS_PATTERN = f"[{BULKY_AAS}]"


# I am assuming here all peptides will have an n-terminal amine and
# a c-terminal carboxy .... aka no, no n- or c- modifications.
def position_index(
    peptide: str,
    pattern: re.Pattern | str = "[KRH]",
    nterm: bool = True,
):
    """Returns the relative position of the first amino acid matching.

    Parameters
    ----------
    peptide : str
        Peptide sequence, no modifications.
    pattern : str | re.Pattern
        Regular expression pattern to match
    nterm : bool
        If True, the position is relative to the n-terminal end of the
        peptide. If False, the position is relative to the c-terminal
        end of the peptide.

    Returns
    -------
    float
        Relative position of the first amino acid matching the pattern.
        If no amino acid matches the pattern, returns 1.
        If it is the first amino acid, returns 0.
        If it is the last amino acid, returns (length(peptide) - 1) / length(peptide).

    Examples
    --------
    >>> position_index("CAAAA", "[C]", nterm=True)
    0.8
    >>> position_index("AAACA", "[C]", nterm=True)
    0.2
    >>> position_index("AAAAAC", "[C]", nterm=True)
    0.0
    >>> position_index("CAAAAA", "[C]", nterm=False)
    0.0
    >>> position_index("LESLIEK", "[KRH]", nterm=False)
    0.8571428571428571
    >>> position_index("LESLIEK", "[KRH]", nterm=True)
    0.0
    """
    if nterm:
        peptide = peptide[::-1]
    m = re.search(pattern, peptide)
    if m:
        return m.start() / len(peptide)
    else:
        return 1


def add_features(
    df,
    stripped_sequence_name="PeptideSequence",
    calc_masses=False,
    default_charge=2,
):
    if np.any(df[stripped_sequence_name].str.contains("[a-z_()0-9]")):
        raise ValueError(
            f"Peptide sequence column {stripped_sequence_name} contains non-amino acid"
            " characters",
        )
    # This two columns are added to later train a simple model that
    # predicts the ion mobility. We will need that model to fill the
    # fields that are missing in .dlib libraries.
    df["StrippedPeptide"] = df[stripped_sequence_name]
    df["PepLength"] = df[stripped_sequence_name].str.len()
    df["NumBulky"] = df[stripped_sequence_name].str.count(BULKY_AAS_PATTERN)
    df["NumPos"] = df[stripped_sequence_name].str.count(POSITIVE_AAS_PATTERN)

    # This gives the relative position of the first [KRH] in the peptide
    # "LESKLIEK" -> 0.25
    # "LESLIEK" -> 0.85
    # "LESLIE" -> 1.0
    df["PosIndexL"] = df[stripped_sequence_name].apply(
        lambda x: position_index(x, POSITIVE_AAS_PATTERN, nterm=False),
    )
    df["PosIndexR"] = df[stripped_sequence_name].apply(
        lambda x: position_index(x, POSITIVE_AAS_PATTERN, nterm=True),
    )

    df["NumNeg"] = df[stripped_sequence_name].str.count(f"[{NEGATIVE_AAS}]")
    df["NegIndexL"] = df[stripped_sequence_name].apply(
        lambda x: position_index(x, NEGATIVE_AAS_PATTERN, nterm=False),
    )
    df["NegIndexR"] = df[stripped_sequence_name].apply(
        lambda x: position_index(x, NEGATIVE_AAS_PATTERN, nterm=True),
    )
    if calc_mass:
        df["Mass"] = df[stripped_sequence_name].apply(
            lambda x: calc_mass(x),
        )
        if "PrecursorCharge" not in df.columns:
            logger.warning(
                "Charge not provided, using default charge of {default_charge}",
            )
            df["PrecursorCharge"] = default_charge

        df["PrecursorMz"] = mass_to_mz(df["Mass"], df["PrecursorCharge"])
    return df


def seq_to_features(stripped_sequence, calc_masses=True, charge=None):
    if any(x in stripped_sequence for x in "[]+_-"):
        raise ValueError(
            "Provided sequence '{stripped_sequence}' contains invalid characters",
        )
    if calc_masses and charge is None:
        raise ValueError("Charge must be provided if calc_masses is True")

    out_features = {}
    out_features["PepLength"] = len(stripped_sequence)
    out_features["NumBulky"] = sum(stripped_sequence.count(x) for x in BULKY_AAS)
    out_features["NumPos"] = sum(stripped_sequence.count(x) for x in POSITIVE_AAS)
    out_features["PosIndexL"] = position_index(
        stripped_sequence,
        POSITIVE_AAS_PATTERN,
        nterm=False,
    )
    out_features["PosIndexR"] = position_index(
        stripped_sequence,
        POSITIVE_AAS_PATTERN,
        nterm=True,
    )
    out_features["NumNeg"] = sum(stripped_sequence.count(x) for x in NEGATIVE_AAS)
    out_features["NegIndexL"] = position_index(
        stripped_sequence,
        NEGATIVE_AAS_PATTERN,
        nterm=False,
    )
    out_features["NegIndexR"] = position_index(
        stripped_sequence,
        NEGATIVE_AAS_PATTERN,
        nterm=True,
    )
    if calc_masses:
        out_features["PrecursorCharge"] = charge
        out_features["Mass"] = calc_mass(stripped_sequence)
        out_features["PrecursorMz"] = mass_to_mz(out_features["Mass"], charge)

    return out_features


def lgb_ims_dataset(df, target_name):
    for col in FEATURE_COLUMNS:
        if col not in df.columns:
            raise ValueError(f"Column {col} not found in dataframe")

    return lgb.Dataset(
        df[FEATURE_COLUMNS],
        label=df[target_name],
        feature_name=FEATURE_COLUMNS,
    )


def df_to_split_datasets(df, target_name, stripped_sequence_name="PeptideSequence"):
    df["Split"] = df[stripped_sequence_name].apply(select_split)
    train_df = df[df["Split"] == "Train"]
    test_df = df[df["Split"] == "Test"]
    val_df = df[df["Split"] == "Val"]

    logger.info(
        f"Split data ({len(df)}) into {len(train_df)} train, {len(val_df)} val,"
        f" {len(test_df)} test",
    )

    train_dataset = lgb_ims_dataset(train_df, target_name)
    val_dataset = lgb_ims_dataset(val_df, target_name)
    test_dataset = lgb_ims_dataset(test_df, target_name)

    return train_dataset, val_dataset, test_dataset
