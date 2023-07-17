# Flimsy: Fun/Fast Simple IMS Anyone like You can use.
Sebastian Paez

This repository implements a very simple LGBM model to predict ion
mobility from peptides.

## Usage

There are two main ways to interact with `flimsay`, one is using python
and the other one is using the python api directly.

### CLI

``` shell
$ pip install flimsay
```

``` shell
$ flimsay fill_blib mylibrary.blib # This will add ion mobility data to a .blib file.
```

``` python
! flimsay fill_blib --help
```


     Usage: flimsay fill_blib [OPTIONS] BLIB OUT_BLIB

     Add ion mobility prediction to a .blib file.

    ╭─ Options ────────────────────────────────────────────────────────────────────╮
    │ --overwrite      Whether to overwrite output file, if it exists              │
    │ --help           Show this message and exit.                                 │
    ╰──────────────────────────────────────────────────────────────────────────────╯

### Python

#### Single peptide

``` python
from flimsay.model import FlimsayModel

model_instance = FlimsayModel()
model_instance.predict_peptide("MYPEPTIDEK", charge=2)
```

    {'ccs': array([363.36245907]), 'one_over_k0': array([0.92423264])}

#### Many peptides at once

``` python
import pandas as pd
from flimsay.features import add_features, FEATURE_COLUMNS

df = pd.DataFrame({
    "Stripped_Seqs": ["LESLIEK", "LESLIE", "LESLKIE"]
})
df = add_features(
    df,
    stripped_sequence_name="Stripped_Seqs",
    calc_masses=True,
    default_charge=2,
)
df
```

    2023-07-17 01:54:10.066 | WARNING  | flimsay.features:add_features:163 - Charge not provided, using default charge of 2

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|     | Stripped_Seqs | StrippedPeptide | PepLength | NumBulky | NumTiny | NumProlines | NumGlycines | NumSerines | NumPos | PosIndexL | PosIndexR | NumNeg | NegIndexL | NegIndexR | Mass       | PrecursorCharge | PrecursorMz |
|-----|---------------|-----------------|-----------|----------|---------|-------------|-------------|------------|--------|-----------|-----------|--------|-----------|-----------|------------|-----------------|-------------|
| 0   | LESLIEK       | LESLIEK         | 7         | 3        | 1       | 0           | 0           | 1          | 1      | 0.857143  | 0.000000  | 2      | 0.142857  | 0.142857  | 830.474934 | 2               | 416.245292  |
| 1   | LESLIE        | LESLIE          | 6         | 3        | 1       | 0           | 0           | 1          | 0      | 1.000000  | 1.000000  | 2      | 0.166667  | 0.000000  | 702.379971 | 2               | 352.197811  |
| 2   | LESLKIE       | LESLKIE         | 7         | 3        | 1       | 0           | 0           | 1          | 1      | 0.571429  | 0.285714  | 2      | 0.142857  | 0.000000  | 830.474934 | 2               | 416.245292  |

</div>

``` python
model_instance.predict(df[FEATURE_COLUMNS])
```

    {'ccs': array([315.32424627, 306.70134752, 314.87268797]),
     'one_over_k0': array([0.78718781, 0.72658194, 0.78525451])}

## Performance

### Prediction Performance

![](https://github.com/TalusBio/flimsay/blob/main/train/plots/one_over_k0_model_ims_pred_vs_true.png)

![](https://github.com/TalusBio/flimsay/blob/main/train/plots/ccs_predicted_vs_real.png)

### Prediction Speed

#### Single peptide prediction

``` python
from flimsay.model import FlimsayModel

model_instance = FlimsayModel()

%timeit model_instance.predict_peptide("MYPEPTIDEK", charge=3)
```

    135 µs ± 10.7 µs per loop (mean ± std. dev. of 7 runs, 10,000 loops each)

In my laptop that takes 133 microseconds per peptide, or roughly 7,500
peptides per second.

#### Batch Prediction

``` python
# Lets make a dataset of 1M peptides to test
import random
import pandas as pd
from flimsay.features import calc_mass, mass_to_mz, add_features

random.seed(42)
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
charges = [2,3,4]

seqs = [random.sample(AMINO_ACIDS, 10) for _ in range(1_000_000)]
charges = [random.sample(charges, 1)[0] for _ in range(1_000_000)]
seqs = ["".join(seq) for seq in seqs]
masses = [calc_mass(x) for x in seqs]
mzs = [mass_to_mz(m, c) for m, c in zip(masses, charges)]

df = pd.DataFrame({
    "Stripped_Seqs": seqs,
    "PrecursorCharge": charges,
    "Mass": masses,
    "PrecursorMz": mzs})
df = add_features(df, stripped_sequence_name="Stripped_Seqs")


# Now we get to run the prediction!
%timeit model_instance.predict(df)
```

    7.79 s ± 1.42 s per loop (mean ± std. dev. of 7 runs, 1 loop each)

In my system every million peptides is predicted in 8.86 seconds, that is
113,000 per second.

## Motivation

There is a fair amount of data on CCS and ion mobility of peptides but
only very few models actually use features that are directly
interpretable.

In addition, having a simpler model allows faster predictions in systems
that are not equiped with GPUs.

Therefore, this project aims to create a freely available, easy to use,
interpretable and fast model to predict ion mobility and collisional
cross-section for peptides.

## Features

The features used for prediction are meant to be simple and their
implementation can be found here:
[flimsy/features.py](flimsy/features.py)

``` python
from flimsay.features import FEATURE_COLUMN_DESCRIPTIONS
for k,v in FEATURE_COLUMN_DESCRIPTIONS.items():
    print(f">>> The Feature '{k}' is defined as: \n\t{v}")
```

    >>> The Feature 'PrecursorMz' is defined as:
        Measured precursor m/z
    >>> The Feature 'Mass' is defined as:
        Measured precursor mass (Da)
    >>> The Feature 'PrecursorCharge' is defined as:
        Measured precursor charge, from the isotope envelope
    >>> The Feature 'PepLength' is defined as:
        Length of the peptide sequence in amino acids
    >>> The Feature 'NumBulky' is defined as:
        Number of bulky amino acids (LVIFWY)
    >>> The Feature 'NumTiny' is defined as:
        Number of tiny amino acids (AS)
    >>> The Feature 'NumProlines' is defined as:
        Number of proline residues
    >>> The Feature 'NumGlycines' is defined as:
        Number of glycine residues
    >>> The Feature 'NumSerines' is defined as:
        Number of serine residues
    >>> The Feature 'NumPos' is defined as:
        Number of positive amino acids (KRH)
    >>> The Feature 'PosIndexL' is defined as:
        Relative position of the first positive amino acid (KRH)
    >>> The Feature 'PosIndexR' is defined as:
        Relative position of the last positive amino acid (KRH)
    >>> The Feature 'NumNeg' is defined as:
        Number of negative amino acids (DE)
    >>> The Feature 'NegIndexL' is defined as:
        Relative position of the first negative amino acid (DE)
    >>> The Feature 'NegIndexR' is defined as:
        Relative position of the last negative amino acid (DE)

## Training

Currently the training logic is handled using DVC (https://dvc.org).

``` shell
git clone {this repo}
cd flimsay/train
dvc repro
```

Running this should automatically download the data, trian the models,
calculate and update the metrics.

The current version of this repo uses predominantly the data from: -
Meier, F., Köhler, N.D., Brunner, AD. et al. Deep learning the
collisional cross sections of the peptide universe from a million
experimental values. Nat Commun 12, 1185 (2021).
https://doi.org/10.1038/s41467-021-21352-8
