from __future__ import annotations

import lightgbm as lgb
import numpy as np
import pandas as pd

from .features import FEATURE_COLUMNS, seq_to_features
from .weights import DEFAULT_CCS_WEIGHTS_PATH, DEFAULT_OOK0_WEIGHTS_PATH


class FlimsayModel:
    """Flimsay model for predicting ion mobility values from peptide sequences."""

    def __init__(self) -> None:
        """Initialize the model using the default weights."""
        self.ccs_model = lgb.Booster(model_file=DEFAULT_CCS_WEIGHTS_PATH)
        self.ook0_model = lgb.Booster(model_file=DEFAULT_OOK0_WEIGHTS_PATH)

    def predict_peptide(self, peptide: str, charge: int = 2) -> dict[str, float]:
        """Predict ion mobility values for a peptide sequence.

        Parameters
        ----------
        peptide : str
            Peptide sequence.
        charge : int, optional
            Peptide charge, by default 2

        Returns
        -------
        dict
            Dictionary of ion mobility values.
            "ccs" and "one_over_k0" are the keys.
        """
        features = seq_to_features(peptide, calc_masses=True, charge=charge)
        feature_array = np.array([[features[x]] for x in FEATURE_COLUMNS]).T

        return self._predict(feature_array)

    def predict(self, df: pd.DataFrame) -> dict[str, np.ndarray | float]:
        """Predict ion mobility values for a dataframe of peptides.

        Parameters
        ----------
        df : pd.DataFrame
            Dataframe of peptides.
            must contain all columns specified in flimsay.features.FEATURE_COLUMNS.
            Check flimsay.features.FEATURE_COLUMNS_DESCRIPTION for more information.
        """
        return self._predict(df[FEATURE_COLUMNS].values)

    def _predict(self, features: np.ndarray) -> dict:
        """Predict ion mobility values for a set of features."""
        ccs = self.ccs_model.predict(features)
        ook0 = self.ook0_model.predict(features)
        return {"ccs": ccs, "one_over_k0": ook0}
