import numpy as np
import pandas as pd
import pytest

from flimsay.features import (
    FEATURE_COLUMN_DESCRIPTIONS,
    add_features,
    calc_mass,
    seq_to_features,
)


def test_add_features():
    """Test that addiiton of features to dataframe works."""
    df = pd.DataFrame({"Stripped_Seqs": ["LESLIEK", "LESLIE", "LESLKIE"]})
    df = add_features(df, stripped_sequence_name="Stripped_Seqs")
    assert df["PepLength"].tolist() == [7, 6, 7]
    assert df["NumBulky"].tolist() == [3, 3, 3]
    assert df["NumPos"].tolist() == [1, 0, 1]
    assert np.allclose(df["PosIndexL"].tolist(), [0.85, 1.0, 0.55], atol=0.1, rtol=0.1)
    assert np.allclose(df["PosIndexR"].tolist(), [0.0, 0.85, 0.25], atol=0.1, rtol=0.1)
    assert df["NumNeg"].tolist() == [2, 2, 2]
    assert np.allclose(
        df["NegIndexL"].tolist(),
        [0.14, 0.16, 0.14],
        atol=0.01,
        rtol=0.01,
    )
    assert np.allclose(df["NegIndexR"].tolist(), [0.14, 0, 0], atol=0.01, rtol=0.01)


def test_seq_to_features():
    """Test that the features are calculated correctly at the peptide level."""
    feats = seq_to_features("LESLIEK")
    assert feats["PepLength"] == 7
    assert feats["NumBulky"] == 3
    assert feats["NumPos"] == 1
    assert abs(feats["PosIndexL"] - 0.85) < 0.01
    assert feats["PosIndexR"] == 0.0
    assert feats["NumNeg"] == 2
    assert abs(feats["NegIndexL"] - 0.14) < 0.01
    assert abs(feats["NegIndexR"] - 0.14) < 0.01


@pytest.fixture
def _sample_sequences():
    """Generate a dataframe of random sequences."""
    import random

    random.seed(42)
    AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")  # noqa
    seqs = [random.sample(AMINO_ACIDS, 10) for _ in range(100)]
    seqs = ["".join(seq) for seq in seqs]
    df = pd.DataFrame({"Stripped_Seqs": seqs})
    return df, seqs


def test_seq_to_features_consistency(_sample_sequences):
    """Test that the features are consistent with the dataframe version."""
    df, seqs = _sample_sequences
    df = add_features(df, stripped_sequence_name="Stripped_Seqs", calc_masses=True)

    for i, seq in enumerate(seqs):
        feats = seq_to_features(seq, calc_masses=True, charge=2)
        for k in FEATURE_COLUMN_DESCRIPTIONS:
            assert df[k].tolist()[i] == feats[k], f"Failed for {k} and {seq}"


def test_masses_fail_with_wrong_aas():
    """Test that calculating masses fail with wrong amino acids."""
    correct_sequences = ["LESLIEK", "LESLIE", "LESLKIE"]
    for seq in correct_sequences:
        calc_mass(seq)

    incorrect_sequences = [
        "_LESLIEK_",
        "LES[UNIMOD:21]LIE",
        "A.LESLKIE.A",
        "ACET-LESLIEK",
        "LES(PHOS)LIEKACET",
    ]
    for seq in incorrect_sequences:
        with pytest.raises(KeyError):
            calc_mass(seq)
