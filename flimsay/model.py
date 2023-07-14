from pathlib import Path

import lightgbm as lgb
import numpy as np

from .features import FEATURE_COLUMNS, seq_to_features

WEIGHTS_PATH = Path(__file__).parent / "weights"
DEFAULT_CCS_WEIGHTS_PATH = WEIGHTS_PATH / "ccs_model.txt"
DEFAULT_OOK0_WEIGHTS_PATH = WEIGHTS_PATH / "one_over_k0_model.txt"


class FlimsayModel:
    def __init__(self) -> None:
        self.ccs_model = lgb.Booster(model_file=DEFAULT_CCS_WEIGHTS_PATH)
        self.ook0_model = lgb.Booster(model_file=DEFAULT_OOK0_WEIGHTS_PATH)

    def predict_peptide(self, peptide: str, charge: int = 2):
        features = seq_to_features(peptide, calc_masses=True, charge=charge)
        feature_array = np.array([[features[x]] for x in FEATURE_COLUMNS]).T

        return self._predict(feature_array)

    def predict(self, df):
        return self._predict(df[FEATURE_COLUMNS].values)

    def _predict(self, features):
        ccs = self.ccs_model.predict(features)
        ook0 = self.ook0_model.predict(features)
        return {"ccs": ccs, "one_over_k0": ook0}
