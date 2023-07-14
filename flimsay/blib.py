import argparse
import sqlite3

import pandas as pd
from loguru import logger

from .features import FEATURE_COLUMNS, add_features, mz_to_mass
from .model import FlimsayModel

pd.set_option("display.max_columns", None)


REF_SPECTRA_SCHEMA = """
CREATE TABLE RefSpectra ( -- spectrum metadata - actual mz/intensity pairs in RefSpectraPeaks
    id INTEGER primary key autoincrement not null, -- lookup key for RefSpectraPeaks
    peptideSeq VARCHAR(150), -- unmodified peptide sequence, can be left blank for small molecule use
    precursorMZ REAL, -- mz of the precursor that produced this spectrum
    precursorCharge INTEGER, -- should agree with adduct if provided
    peptideModSeq VARCHAR(200), -- modified peptide sequence, can be left blank for small molecule use
    prevAA CHAR(1), -- position of peptide in its parent protein (can be left blank)
    nextAA CHAR(1),  -- position of peptide in its parent protein (can be left blank)
    copies INTEGER, -- number of copies this spectrum was chosen from if it is in a filtered library
    numPeaks INTEGER, -- number of peaks, should agree with corresponding entry in RefSpectraPeaks
    ionMobility REAL, -- ion mobility value, if known (see ionMobilityType for units)
    collisionalCrossSectionSqA REAL, -- precursor CCS in square Angstroms for ion mobility, if known
    ionMobilityHighEnergyOffset REAL, -- ion mobility value increment for fragments (see ionMobilityType for units)
    ionMobilityType TINYINT, -- ion mobility units (required if ionMobility is used, see IonMobilityTypes table for key)
    retentionTime REAL, -- chromatographic retention time in minutes, if known
    startTime REAL, -- start retention time in minutes, if known
    endTime REAL, -- end retention time in minutes, if known
    totalIonCurrent REAL, -- total ion current of spectrum
    moleculeName VARCHAR(128), -- precursor molecule's name (not needed for peptides)
    chemicalFormula VARCHAR(128), -- precursor molecule's neutral formula (not needed for peptides)
    precursorAdduct VARCHAR(128), -- ionizing adduct e.g. [M+Na], [2M-H2O+2H] etc (not needed for peptides)
    inchiKey VARCHAR(128), -- molecular identifier for structure retrieval (not needed for peptides)
    otherKeys VARCHAR(128), -- alternative molecular identifiers for structure retrieval, tab separated name:value pairs e.g. cas:58-08-2\thmdb:01847 (not needed for peptides)
    fileID INTEGER, -- index into SpectrumSourceFiles table for source file information
    SpecIDinFile VARCHAR(256), -- original spectrum label, id, or description in source file
    score REAL, -- spectrum score, typically a probability score (see scoreType)
    scoreType TINYINT -- spectrum score type, see ScoreTypes table for meaning
);
CREATE INDEX idxPeptide ON RefSpectra (peptideSeq, precursorCharge);
CREATE INDEX idxPeptideMod ON RefSpectra (peptideModSeq, precursorCharge);
CREATE INDEX idxMoleculeName ON RefSpectra (moleculeName, precursorAdduct);
CREATE INDEX idxInChiKey ON RefSpectra (inchiKey, precursorAdduct);
"""  # noqa

RT_SCHEMA = """
CREATE TABLE IF NOT EXISTS "RetentionTimes" (
    "RefSpectraID" INTEGER,
    "RedundantRefSpectraID" INTEGER,
    "SpectrumSourceID" INTEGER,
    "ionMobility" REAL,
    "collisionalCrossSectionSqA" REAL,
    "ionMobilityHighEnergyOffset" REAL,
    "ionMobilityType" INTEGER,
    "retentionTime" REAL,
    "startTime" TEXT,
    "endTime" TEXT,
    "score" REAL,
    "bestSpectrum" INTEGER
);
"""

LIBINFO_SCHEMA = """
CREATE TABLE LibInfo(
    libLSID TEXT,
    createTime TEXT,
    numSpecs INTEGER,
    majorVersion INTEGER,
    minorVersion INTEGER
);
"""


def main(blib):
    pred_model = FlimsayModel()
    conn = sqlite3.connect(blib)
    cur = conn.cursor()
    logger.info("Reading from file")
    df = pd.read_sql_query("SELECT * FROM RefSpectra", conn)
    logger.info(df.head())

    pred_df = df[["id", "peptideSeq", "precursorMZ", "precursorCharge"]].copy()
    pred_df["PeptideSequence"] = pred_df["peptideSeq"]
    pred_df["PrecursorMz"] = pred_df["precursorMZ"]
    pred_df["PrecursorCharge"] = pred_df["precursorCharge"]
    pred_df["Mass"] = mz_to_mass(pred_df["PrecursorMz"], pred_df["PrecursorCharge"])

    pred_df = add_features(pred_df, calc_masses=False)
    missing_cols = set(FEATURE_COLUMNS) - set(pred_df.columns)
    if missing_cols:
        logger.error("Missing columns: ", missing_cols)
        logger.error("Columns: ", pred_df.columns)
        raise KeyError("Missing columns")

    logger.info("Predicting ion mobility")
    ######### Prediction Start
    predictions = pred_model.predict(pred_df)
    ook0_predicted = predictions["one_over_k0"]
    ccs_predicted = predictions["ccs"]

    id_to_imns = dict(zip(pred_df["id"], ook0_predicted))
    id_to_ccs = dict(zip(pred_df["id"], ccs_predicted))
    ####### Prediction End

    try:
        del df["ionMobilityValue"]
    except KeyError:
        pass
    df["collisionalCrossSectionSqA"] = ccs_predicted
    df["ionMobility"] = ook0_predicted
    df["ionMobilityHighEnergyOffset"] = 0
    im_offset = df["ionMobilityHighEnergyOffset"].astype("float64")
    df["ionMobilityHighEnergyOffset"] = im_offset
    df["ionMobilityType"] = 2

    logger.info("Writing to file")
    logger.info("Updating RefSpectra Table")

    cur.execute("DROP TABLE RefSpectra")
    cur.executescript(REF_SPECTRA_SCHEMA)
    df.to_sql(
        "RefSpectra",
        conn,
        if_exists="append",
        index=False,
        schema=REF_SPECTRA_SCHEMA,
    )

    logger.info("Updating RetentionTimes Table")
    rtdf = pd.read_sql_query("SELECT * FROM RetentionTimes", conn)

    try:
        del rtdf["ionMobilityValue"]
    except KeyError:
        pass

    rtdf["ionMobility"] = [id_to_imns[x] for x in rtdf["RefSpectraID"]]
    rtdf["collisionalCrossSectionSqA"] = [id_to_ccs[x] for x in rtdf["RefSpectraID"]]
    rtdf["ionMobilityType"] = 2
    rtdf["ionMobilityHighEnergyOffset"] = float(0)

    # Making sure this is a float
    im_offset = rtdf["ionMobilityHighEnergyOffset"].astype("float64")
    rtdf["ionMobilityHighEnergyOffset"] = im_offset

    cur.execute("DROP TABLE RetentionTimes")
    cur.executescript(RT_SCHEMA)
    rtdf.to_sql(
        "RetentionTimes",
        conn,
        if_exists="replace",
        index=False,
        schema=RT_SCHEMA,
    )

    # Just a print for sanity checking
    logger.debug("Printing header of the new table")
    df = pd.read_sql_query("SELECT * FROM RefSpectra LIMIT 5", conn)
    logger.debug(df.head())

    df = pd.read_sql_query("SELECT * FROM RetentionTimes LIMIT 5", conn)
    logger.debug(df.head())

    logger.info("Creating ims info table")
    # Encyclopedia does not write this table, so we add it here
    try:
        cur.executescript("""
            CREATE TABLE IonMobilityTypes (
                id INTEGER PRIMARY KEY,
                ionMobilityType VARCHAR(128)
            );
            INSERT INTO IonMobilityTypes(id, ionMobilityType)
                VALUES(0, 'none');
            INSERT INTO IonMobilityTypes(id, ionMobilityType)
                VALUES(1, 'driftTime(msec)');
            INSERT INTO IonMobilityTypes(id, ionMobilityType)
                VALUES(2, 'inverseK0(Vsec/cm^2)');
            INSERT INTO IonMobilityTypes(id, ionMobilityType)
                VALUES(3, 'compensation(V)');
            """)
    except sqlite3.OperationalError as e:
        if "table IonMobilityTypes already exists" in str(e):
            pass
        else:
            raise

    logger.info("Modifying lib info table")
    df = pd.read_sql_query("SELECT * FROM LibInfo", conn)
    df["minorVersion"] = 6

    cur.execute("DROP TABLE LibInfo")
    cur.executescript(LIBINFO_SCHEMA)
    df.to_sql("LibInfo", conn, if_exists="append", index=False, schema=LIBINFO_SCHEMA)

    conn.commit()
    conn.close()
    logger.info("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add IMS to BLIB")
    parser.add_argument("blib", help="BLIB file")

    args, unknown = parser.parse_known_args()
    if unknown:
        raise RuntimeError("Unrecognized arguments: ", unknown)
    else:
        main(
            args.blib,
            one_over_k0_model_file=args.one_over_k0_model_file,
            ccs_model_file=args.ccs_model_file,
        )
