import re

import lightgbm as lgb
import numpy as np
import pandas as pd
from loguru import logger

from .data import select_split

FEATURE_COLUMNS = [
    "PrecursorMz",
    "Mass",
    "PrecursorCharge",
    "PepLength",
    "NumBulky",
    "NumPos",
    "PosIndexL",
    "PosIndexR",
    "NumNeg",
    "NegIndexL",
    "NegIndexR",
]

BULKY_AAS = "LVIFWY"
POSITIVE_AAS = "KRH"
NEGATIVE_AAS = "DE"

POSITIVE_AAS_PATTERN = f"[{POSITIVE_AAS}]"
NEGATIVE_AAS_PATTERN = f"[{NEGATIVE_AAS}]"
BULKY_AAS_PATTERN = f"[{BULKY_AAS}]"


# I am assuming here all peptides will have an n-terminal amine and
# a c-terminal carboxy .... aka no, no n- or c- modifications.
def position_index(peptide, pattern="[KRH]", nterm=True):
    if nterm:
        peptide = peptide[::-1]
    m = re.search(pattern, peptide)
    if m:
        return m.start() / len(peptide)
    else:
        return 1


def add_features(df, stripped_sequence_name="PeptideSequence"):
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
    return df


def test_add_features():
    df = pd.DataFrame({"Stripped_Seqs": ["LESLIEK", "LESLIE", "LESLKIE"]})
    df = add_features(df, stripped_sequence_name="Stripped_Seqs")
    assert df["PepLength"].tolist() == [7, 6, 7]
    assert df["NumBulky"].tolist() == [3,3,3]
    assert df["NumPos"].tolist() == [1, 0, 1]
    assert np.allclose(df["PosIndexL"].tolist(), [0.85, 1.0, 0.55], atol=0.1, rtol=0.1)
    assert np.allclose(df["PosIndexR"].tolist(), [0.0, 0.85, 0.25], atol=0.1, rtol=0.1)
    assert df["NumNeg"].tolist() == [2, 2, 2]
    assert np.allclose(df["NegIndexL"].tolist(), [0.14, 0.16, 0.14], atol=0.01, rtol=0.01)
    assert np.allclose(df["NegIndexR"].tolist(), [0.14, 0, 0], atol=0.01, rtol=0.01)


def seq_to_features(stripped_sequence):
    if any(x in stripped_sequence for x in "[]+_-"):
        raise ValueError("Provided sequence '{stripped_sequence}' contains invalid characters")
    out_features = {}
    out_features["PepLength"] = len(stripped_sequence)
    out_features["NumBulky"] = sum(stripped_sequence.count(x) for x in BULKY_AAS)
    out_features["NumPos"] = sum(stripped_sequence.count(x) for x in POSITIVE_AAS)
    out_features["PosIndexL"] = position_index(
        stripped_sequence, POSITIVE_AAS_PATTERN, nterm=False,
    )
    out_features["PosIndexR"] = position_index(
        stripped_sequence, POSITIVE_AAS_PATTERN, nterm=True,
    )
    out_features["NumNeg"] = sum(stripped_sequence.count(x) for x in NEGATIVE_AAS)
    out_features["NegIndexL"] = position_index(
        stripped_sequence, NEGATIVE_AAS_PATTERN, nterm=False,
    )
    out_features["NegIndexR"] = position_index(
        stripped_sequence, NEGATIVE_AAS_PATTERN, nterm=True,
    )
    return out_features


def test_seq_to_features():
    feats = seq_to_features("LESLIEK")
    assert feats["PepLength"] == 7
    assert feats["NumBulky"] == 3
    assert feats["NumPos"] == 1
    assert abs(feats["PosIndexL"] - 0.85) < 0.01
    assert feats["PosIndexR"] == 0.0
    assert feats["NumNeg"] == 2
    assert abs(feats["NegIndexL"] - 0.14) < 0.01
    assert abs(feats["NegIndexR"] - 0.14) < 0.01


def test_seq_to_features_consistency():
    import random

    random.seed(42)
    AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
    seqs = [random.sample(AMINO_ACIDS, 10) for _ in range(100)]
    seqs = ["".join(seq) for seq in seqs]
    df = pd.DataFrame({"Stripped_Seqs": seqs})
    df = add_features(df, stripped_sequence_name="Stripped_Seqs")

    for i, seq in enumerate(seqs):
        feats = seq_to_features(seq)
        for k, v in feats.items():
            assert df[k].tolist()[i] == v, f"Failed for {k} and {seq}"


def mz_to_mass(mz, charge):
    return (mz * charge) - charge * 1.00727646677


def mass_to_mz(mass, charge):
    return (mass + charge * 1.00727646677) / charge


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
