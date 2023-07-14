import numpy as np
import pandas as pd

from flimsay.features import add_features, seq_to_features


def test_add_features():
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
