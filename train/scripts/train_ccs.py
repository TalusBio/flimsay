from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
import vizta
from dataset_utils import find_mod_pairs, write_simple_toml
from lightgbm.callback import log_evaluation
from matplotlib import pyplot as plt
from matplotlib import rcParams

from flimsay.data import select_split
from flimsay.features import (
    FEATURE_COLUMNS,
    add_features,
    df_to_split_datasets,
    mass_to_mz,
)

WEIGHTS_LOCATION = Path(__file__).parent / "../weights"
PLOTS_LOCATION = Path(__file__).parent / "../plots"
DATA_LOCATION = Path(__file__).parent / "../data"

model_name = "ccs_model.txt"
model_location = WEIGHTS_LOCATION / model_name

vizta.mpl.set_theme()  # Set the theme
rcParams.update({"figure.autolayout": True})

peplib = pd.read_csv(DATA_LOCATION / "ccs_data.csv", index_col=0)

peplib["Modified sequence"] = peplib["Modified sequence"].str.replace("_", "")
peplib["PeptideSequence"] = peplib["Modified sequence"].str.replace(
    "[a-z_()]",
    "",
    regex=True,
)


_ = find_mod_pairs(
    peplib["Modified sequence"].to_list(),
    peplib["PeptideSequence"].to_list(),
)
peplib["PrecursorCharge"] = peplib["Charge"]
peplib["PrecursorMz"] = mass_to_mz(peplib["Mass"], peplib["Charge"])

peplib = add_features(peplib, calc_masses=False)
train_data, validation_data, test_data = df_to_split_datasets(
    peplib,
    "CCS",
    "PeptideSequence",
)
peplib["split"] = peplib["PeptideSequence"].apply(select_split)


num_round = 2000
param = {
    "objective": "regression",
    "metric": ["l1", "l2_root"],
    "num_leaves": 51,
    "early_stopping_rounds": 20,
    "learning_rate": 0.1,
    "pred_early_stop_margin": 1e-3,
}

bst = lgb.train(
    param,
    train_set=train_data,
    num_boost_round=num_round,
    valid_sets=[validation_data],
    callbacks=[log_evaluation(period=20)],
)

bst.save_model(model_location, num_iteration=bst.best_iteration)
validation_data = peplib[peplib["split"] == "Val"].copy()
validation_data = add_features(validation_data)

predicted = bst.predict(validation_data[FEATURE_COLUMNS])
validation_data["predicted"] = predicted
validation_data["error"] = validation_data["CCS"] - predicted
metrics = {}
avg_error = validation_data["error"].abs().mean()
metrics["CCS_MAE"] = avg_error

plt.scatter(
    validation_data["PrecursorMz"],
    validation_data["CCS"] - predicted,
    s=1,
    c=validation_data["PrecursorCharge"],
    alpha=0.2,
)
# plt.ylim(-0.10, 0.10)
plt.xlabel("PrecursorMz")
plt.ylabel("CCS - predicted CCS")
plt.colorbar()
plt.savefig(PLOTS_LOCATION / "ccs_mz_vs_error.png")
plt.close()

plt.hist(validation_data["error"], bins=100)
plt.xlabel("CCS - predicted CCS")
plt.ylabel("Count")
plt.title(
    "Validation set CCS prediction error"
    f" MAE={validation_data['error'].abs().mean():.04f}\n{model_name}",
)
plt.savefig(PLOTS_LOCATION / "ccs_error.png")
plt.close()


m, b = np.polyfit(validation_data["CCS"], predicted, deg=1)
pearson = np.corrcoef(validation_data["CCS"], predicted, 1)[0, 1]
limits = [validation_data["CCS"].min(), validation_data["CCS"].max()]
plt.axline(xy1=(0, b), slope=m, label=f"$y = {m:.1f}x {b:+.1f}$", c="#666666")
metrics["CCS_Pearson"] = pearson

plt.scatter(
    predicted,
    validation_data["CCS"],
    s=1,
    c=validation_data["PrecursorCharge"],
    alpha=0.2,
)
plt.xlabel("Predicted CCS")
plt.ylabel("CCS")
plt.title(f"Predicted vs Real CCS (Pearson: {float(pearson):.4f})")
plt.xlim(*limits)
plt.ylim(*limits)
plt.colorbar()
plt.savefig(PLOTS_LOCATION / "ccs_predicted_vs_real.png")
plt.close()

lgb.plot_importance(bst, figsize=(6, 3))
plt.savefig(PLOTS_LOCATION / "ccs_feature_importance.png")
plt.close()

# Plot with y axis being the cumulative number of peptides
# and y the error
plt.plot(
    np.sort(np.abs(validation_data["error"])),
    np.arange(0, len(validation_data)),
    label="All",
)
# label what percentage is under 20 error
num_under_20 = len(validation_data[validation_data["error"] < 20])
pct_under_20 = num_under_20 / len(validation_data)
plt.text(
    20,
    num_under_20,
    f"{pct_under_20:.2%} under 20 error",
    horizontalalignment="center",
)
for i in np.unique(validation_data["PrecursorCharge"]):
    plt.plot(
        np.sort(
            np.abs(validation_data[validation_data["PrecursorCharge"] == i]["error"]),
        ),
        np.arange(0, len(validation_data[validation_data["PrecursorCharge"] == i])),
        label=f"Charge {i}",
    )
    # label what percentage of the data is within 0.01
    # of the predicted value
    num_under_20 = len(
        validation_data[validation_data["PrecursorCharge"] == i][
            np.abs(validation_data[validation_data["PrecursorCharge"] == i]["error"])
            < 20
        ],
    )
    pct_under_20 = num_under_20 / len(
        validation_data[validation_data["PrecursorCharge"] == i],
    )
    plt.text(
        20,
        num_under_20,
        f"{pct_under_20:.2f}",
        horizontalalignment="left",
        verticalalignment="center",
    )

# Show only until 100
plt.xlim(0, 50)
plt.legend()
plt.ylabel("Cumulative Count")
plt.xlabel("abs(IonMobility - predicted Mobility)")
plt.title(f"Validation Data MAE: {avg_error:.4f}")
plt.savefig(PLOTS_LOCATION / "ccs_cumulative_error.png")

# Save metrics
# This is a pretty simple and plain toml file
write_simple_toml(PLOTS_LOCATION / "ccs_model_metrics.toml", metrics)
