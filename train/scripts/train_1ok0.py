from pathlib import Path

import lightgbm as lgb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vizta
from dataset_utils import find_mod_pairs, write_simple_toml
from lightgbm.callback import log_evaluation
from loguru import logger
from matplotlib import rcParams

from flimsay.features import (
    FEATURE_COLUMNS,
    add_features,
    df_to_split_datasets,
    mz_to_mass,
    select_split,
)

vizta.mpl.set_theme()  # Set the theme
rcParams.update({"figure.autolayout": True})

## Columns used from the data
# PrecursorCharge <<
# ModifiedPeptide <<
# StrippedPeptide <<
# iRT <<
# IonMobility <<
# PrecursorMz <<

# All columsn in the data
# ReferenceRun
# PrecursorCharge <<
# Workflow
# IntModifiedPeptide
# CV
# AllowForNormalization
# ModifiedPeptide <<
# StrippedPeptide <<
# iRT <<
# IonMobility <<
# iRTSourceSpecific
# BGSInferenceId
# IsProteotypic
# IntLabeledPeptide
# LabeledPeptide
# PrecursorMz <<
# ReferenceRunQvalue
# ReferenceRunMS1Response
# FragmentLossType
# FragmentNumber
# FragmentType
# FragmentCharge
# FragmentMz
# RelativeIntensity
# ExcludeFromAssay
# Database
# ProteinGroups
# UniProtIds
# Protein
# Name
# ProteinDescription
# Organisms
# OrganismId
# Genes
# Protein
# Existence
# Sequence
# Version
# FASTAName

WEIGHTS_LOCATION = Path(__file__).parent / "../weights"
PLOTS_LOCATION = Path(__file__).parent / "../plots"

peplib = pd.read_csv(
    # Path(__file__).parent / "../IM-GPF_library.tsv",
    Path(__file__).parent / "../data/Deep_Library.tsv",
    sep="\t",
    usecols=[
        "ModifiedPeptide",
        "StrippedPeptide",
        "PrecursorCharge",
        "PrecursorMz",
        "iRT",
        "IonMobility",
    ],
).drop_duplicates(keep="first")

peplib["ModifiedPeptide"] = peplib["ModifiedPeptide"].str.replace("_", "")
peplib = add_features(peplib, stripped_sequence_name="StrippedPeptide")
peplib["Mass"] = mz_to_mass(peplib["PrecursorMz"], peplib["PrecursorCharge"])
logger.info(peplib)


_ = find_mod_pairs(
    peplib["ModifiedPeptide"].to_list(),
    peplib["StrippedPeptide"].to_list(),
)

train_dataset, val_dataset, test_dataset = df_to_split_datasets(
    peplib,
    "IonMobility",
    "StrippedPeptide",
)

model_name = "one_over_k0_model.txt"
model_location = WEIGHTS_LOCATION / model_name

num_round = 2000
param = {
    "objective": "regression",
    "metric": ["l1", "l2_root"],
    "num_leaves": 51,
    "early_stopping_rounds": 20,
    "learning_rate": 0.05,
    "pred_early_stop_margin": 1e-4,
}

bst = lgb.train(
    param,
    train_dataset,
    num_round,
    valid_sets=[val_dataset],
    callbacks=[log_evaluation(period=20)],
)
bst.save_model(model_location, num_iteration=bst.best_iteration)

lgb.plot_importance(bst, figsize=(6, 3))
plt.savefig(PLOTS_LOCATION / "one_over_k0_model_feature_importance.png")
plt.close()


def gen_eval_plots(use_df, ds_name="test"):
    metrics = {}
    use_df["Predicted"] = bst.predict(use_df[FEATURE_COLUMNS])
    plt.scatter(
        use_df["PrecursorMz"],
        use_df["IonMobility"] - use_df["Predicted"],
        s=1,
        c=use_df["PrecursorCharge"],
        alpha=0.1,
    )
    plt.ylim(-0.10, 0.10)
    plt.xlabel("PrecursorMz")
    plt.ylabel("IonMobility - predicted Mobility")
    plt.savefig(PLOTS_LOCATION / "one_over_k0_model_ims_error.png")
    plt.close()

    use_df["IonMobilityError"] = use_df["IonMobility"] - use_df["Predicted"]
    avg_error = np.mean(np.abs(use_df["IonMobilityError"]))
    metrics["OOK0_MAE"] = avg_error
    bins = np.linspace(
        use_df["IonMobilityError"].min(),
        use_df["IonMobilityError"].max(),
        100,
    )
    plt.hist(
        use_df["IonMobility"] - use_df["Predicted"],
        bins=bins,
        label="All",
        alpha=0.2,
    )
    for i in np.unique(use_df["PrecursorCharge"]):
        plt.hist(
            use_df[use_df["PrecursorCharge"] == i]["IonMobilityError"],
            bins=bins,
            label=f"Charge {i}",
            alpha=0.5,
        )

    plt.legend()
    plt.ylabel("Count")
    plt.xlabel("IonMobility - predicted Mobility")
    plt.title(f"{ds_name} MAE: {avg_error:.4f}")
    plt.yscale("log")
    plt.savefig(PLOTS_LOCATION / "one_over_k0_model_ims_error_hist.png")
    plt.close()

    # Plot with y axis being the cumulative number of peptides
    # and y the error
    plt.plot(
        np.sort(np.abs(use_df["IonMobilityError"])),
        np.arange(0, len(use_df)),
        label="All",
    )
    # label what percentage is under 3% error
    num_under_003 = len(use_df[use_df["IonMobilityError"] < 0.03])
    pct_under_003 = num_under_003 / len(use_df)
    plt.text(
        0.03,
        num_under_003,
        f"{pct_under_003:.2%} under 3% error",
        horizontalalignment="center",
    )
    for i in np.unique(use_df["PrecursorCharge"]):
        plt.plot(
            np.sort(np.abs(use_df[use_df["PrecursorCharge"] == i]["IonMobilityError"])),
            np.arange(0, len(use_df[use_df["PrecursorCharge"] == i])),
            label=f"Charge {i}",
        )
        # label what percentage of the data is within 0.01
        # of the predicted value
        num_under_001 = len(
            use_df[use_df["PrecursorCharge"] == i][
                np.abs(use_df[use_df["PrecursorCharge"] == i]["IonMobilityError"])
                < 0.01
            ],
        )
        pct_under_001 = num_under_001 / len(use_df[use_df["PrecursorCharge"] == i])
        plt.text(
            0.01,
            num_under_001,
            f"{pct_under_001:.2f}",
            horizontalalignment="left",
            verticalalignment="center",
        )

    # Show only until 0.1
    plt.xlim(0, 0.1)
    plt.legend()
    plt.ylabel("Cumulative Count")
    plt.xlabel("abs(IonMobility - predicted Mobility)")
    plt.title(f"{ds_name} MAE: {avg_error:.4f}")
    plt.savefig(PLOTS_LOCATION / "one_over_k0_model_ims_error_cumulative.png")
    plt.close()

    m, b = np.polyfit(use_df["IonMobility"], use_df["Predicted"], deg=1)
    pearson = np.corrcoef(use_df["IonMobility"], use_df["Predicted"])[0, 1]
    metrics["OOK0_Pearson"] = pearson
    plt.scatter(
        use_df["Predicted"],
        use_df["IonMobility"],
        s=1,
        c=use_df["PrecursorCharge"],
        alpha=0.2,
    )
    plt.axline(xy1=(0, b), slope=m, label=f"$y = {m:.1f}x {b:+.1f}$", c="#666666")
    plt.xlabel("Predicted Mobility")
    plt.ylabel("IonMobility")
    plt.colorbar(alpha=1)
    plt.xlim(use_df["Predicted"].min(), use_df["Predicted"].max())
    plt.ylim(use_df["IonMobility"].min(), use_df["IonMobility"].max())
    plt.title(f"Predicted vs IonMobility (Pearson: {pearson:.4f})\n{ds_name}")
    plt.savefig(PLOTS_LOCATION / "one_over_k0_model_ims_pred_vs_true.png")
    plt.close()

    return metrics


peplib["Split"] = peplib["StrippedPeptide"].apply(select_split)
val_df = peplib[peplib["Split"] == "Val"].copy()
metrics = gen_eval_plots(val_df, "Validation DS")
logger.info(metrics)


# Save metrics
write_simple_toml(PLOTS_LOCATION / "one_over_k0_model_metrics.toml", metrics)
