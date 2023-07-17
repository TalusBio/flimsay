from io import StringIO

import numpy as np
import pandas as pd
import pytest
from loguru import logger

from flimsay.features import add_features, mz_to_mass
from flimsay.model import FlimsayModel

GT_DATA = """
746635	3	_EVGEHLVSIKK_	EVGEHLVSIKK	-8.742136	0.717	413.574955	Filamin-B
556395	4	_YAPTEVGLHEMHIK_	YAPTEVGLHEMHIK	21.931534	0.763	406.959255	Filamin-B
660946	3	_FADEHVPGSPFTVK_	FADEHVPGSPFTVK	35.291040	0.771	510.924424	Filamin-B
660289	3	_LPNNHIGISFIPR_	LPNNHIGISFIPR	46.744590	0.771	493.280654	Filamin-B
392153	2	_AYGPGLEK_	AYGPGLEK	-14.577010	0.781	417.721436	Filamin-B
583871	1	_SPFEVSVDK_	SPFEVSVDK	18.751308	1.264	1007.504403	Filamin-B
460900	2	_DAGEGLLAVQITDQEGKPK_	DAGEGLLAVQITDQEGKPK	48.206657	1.285	985.015469	Filamin-B
667983	2	_DGTC[Carbamidomethyl (C)]TVTYLPTLPGDYSILVK_	DGTCTVTYLPTLPGDYSILVK	128.502300	1.292	1157.087777	Filamin-B
25821	2	_DGTYAVTYVPLTAGMYTLTMK_	DGTYAVTYVPLTAGMYTLTMK	146.747960	1.298	1148.565623	Filamin-B
562396	2	_DGTYAVTYVPLTAGMYTLTM[Oxidation (M)]K_	DGTYAVTYVPLTAGMYTLTMK	96.014000	1.305	1156.563080	Filamin-B
"""  # noqa

COLNAMES = "	PrecursorCharge	ModifiedPeptide	StrippedPeptide	iRT	IonMobility	PrecursorMz	ProteinDescription"  # noqa


@pytest.fixture
def _gt_data():
    return pd.read_csv(StringIO(GT_DATA), sep="\t", header=None, names=COLNAMES.split())


def test_predictions_stripped_peptide(_gt_data):
    """Tests that predictions are close to ground truth using the model weights."""
    model = FlimsayModel()
    tmp_iter = zip(
        _gt_data["StrippedPeptide"].to_list(),
        _gt_data["IonMobility"].to_list(),
        _gt_data["PrecursorCharge"].to_list(),
    )
    for seq, gt_ims, charge in tmp_iter:
        pred = model.predict_peptide(seq, charge)
        assert pred["one_over_k0"] == pytest.approx(gt_ims, abs=0.08)


def test_predictions_batched(_gt_data):
    """Tests that predictions are close to ground truth using the model weights."""
    model = FlimsayModel()
    df = add_features(_gt_data, stripped_sequence_name="StrippedPeptide")
    df["Mass"] = mz_to_mass(df["PrecursorMz"], df["PrecursorCharge"])
    preds = model.predict(df)

    closeness = np.allclose(preds["one_over_k0"], _gt_data["IonMobility"], atol=0.08)
    mae = np.mean(np.abs(preds["one_over_k0"] - _gt_data["IonMobility"]))
    logger.info(f"MAE: {mae}")
    assert closeness
