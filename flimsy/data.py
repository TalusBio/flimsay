from typing import Literal

IRT_PEPTIDES = {
    "LGGNEQVTR": {"vendor": "biognosys", "irt": -24.92},
    "GAGSSEPVTGLDAK": {"vendor": "biognosys", "irt": 0},
    "VEATFGVDESNAK": {"vendor": "biognosys", "irt": 12.39},
    "YILAGVENSK": {"vendor": "biognosys", "irt": 19.79},
    "TPVISGGPYEYR": {"vendor": "biognosys", "irt": 28.71},
    "TPVITGAPYEYR": {"vendor": "biognosys", "irt": 33.38},
    "DGLDAASYYAPVR": {"vendor": "biognosys", "irt": 42.26},
    "ADVTPADFSEWSK": {"vendor": "biognosys", "irt": 54.62},
    "GTFIIDPGGVIR": {"vendor": "biognosys", "irt": 70.52},
    "GTFIIDPAAVIR": {"vendor": "biognosys", "irt": 87.23},
    "LFLQFGAQGSPFLK": {"vendor": "biognosys", "irt": 100},
    "HEHISSDYAGK": {"vendor": "procal", "irt": -36.83},
    "IGYDHGHIEHK": {"vendor": "procal", "irt": -33.5},
    "TFAHTESHISK": {"vendor": "procal", "irt": -33.32},
    "ISLGEHEGGGK": {"vendor": "procal", "irt": -18.54},
    "YVGDSYDSSAK": {"vendor": "procal", "irt": -16.87},
    "FGTGTYAGGEK": {"vendor": "procal", "irt": -9.35},
    "LSSGYDGTSYK": {"vendor": "procal", "irt": -8.82},
    "TASGVGGFSTK": {"vendor": "procal", "irt": -4.18},
    "LTSGDFGEDSK": {"vendor": "procal", "irt": -3.76},
    "AGDEALGDTYK": {"vendor": "procal", "irt": -3.52},
    "SYASDFGSSAK": {"vendor": "procal", "irt": 1.79},
    "LYSYYSSTESK": {"vendor": "procal", "irt": 6.39},
    "FASDTSDEAFK": {"vendor": "procal", "irt": 7.2},
    "LTDTFADDDTK": {"vendor": "procal", "irt": 8.25},
    "LYTGAGYDEVK": {"vendor": "procal", "irt": 10.53},
    "TLIAYDDSSTK": {"vendor": "procal", "irt": 14.98},
    "TASEFDSAIAQDK": {"vendor": "procal", "irt": 17.84},
    "HDLDYGIDSYK": {"vendor": "procal", "irt": 19.86},
    "FLASSEGGFTK": {"vendor": "procal", "irt": 20.88},
    "HTAYSDFLSDK": {"vendor": "procal", "irt": 25.9},
    "FVGTEYDGLAK": {"vendor": "procal", "irt": 26.82},
    "YALDSYSLSSK": {"vendor": "procal", "irt": 32},
    "YYGTIEDTEFK": {"vendor": "procal", "irt": 33.73},
    "GFLDYESTGAK": {"vendor": "procal", "irt": 35.9},
    "HLTGLTFDTYK": {"vendor": "procal", "irt": 36.5},
    "YFGYTSDTFGK": {"vendor": "procal", "irt": 41.42},
    "HDTVFGSYLYK": {"vendor": "procal", "irt": 41.42},
    "FSYDGFEEDYK": {"vendor": "procal", "irt": 44.22},
    "ALFSSITDSEK": {"vendor": "procal", "irt": 44.88},
    "LYLSEYDTIGK": {"vendor": "procal", "irt": 48.16},
    "HFALFSTDVTK": {"vendor": "procal", "irt": 50.41},
    "VSGFSDISIYK": {"vendor": "procal", "irt": 51.67},
    "GSGGFTEFDLK": {"vendor": "procal", "irt": 51.97},
    "TFTGTTDSFFK": {"vendor": "procal", "irt": 52.2},
    "TFGTETFDTFK": {"vendor": "procal", "irt": 54.53},
    "YTSFYGAYFEK": {"vendor": "procal", "irt": 56.65},
    "LTDELLSEYYK": {"vendor": "procal", "irt": 57.66},
    "ASDLLSGYYIK": {"vendor": "procal", "irt": 57.68},
    "YGFSSEDIFTK": {"vendor": "procal", "irt": 57.77},
    "HTYDDEFFTFK": {"vendor": "procal", "irt": 58.44},
    "FLFTGYDTSVK": {"vendor": "procal", "irt": 61.07},
    "GLSDYLVSTVK": {"vendor": "procal", "irt": 61.34},
    "VYAETLSGFIK": {"vendor": "procal", "irt": 62.57},
    "GLFYGGYEFTK": {"vendor": "procal", "irt": 62.96},
    "GSTDDGFIILK": {"vendor": "procal", "irt": 63.07},
    "TSIDSFIDSYK": {"vendor": "procal", "irt": 63.51},
    "TLLLDAEGFEK": {"vendor": "procal", "irt": 65.49},
    "GFVIDDGLITK": {"vendor": "procal", "irt": 66.46},
    "GFEYSIDYFSK": {"vendor": "procal", "irt": 66.9},
    "GIFGAFTDDYK": {"vendor": "procal", "irt": 71.49},
    "LEIYTDFDAIK": {"vendor": "procal", "irt": 71.99},
    "FTEGGILDLYK": {"vendor": "procal", "irt": 72.95},
    "LLFSYSSGFVK": {"vendor": "procal", "irt": 73.23},
    "STFFSFGDVGK": {"vendor": "procal", "irt": 74.29},
    "LTAYFEDLELK": {"vendor": "procal", "irt": 75.09},
    "VDTFLDGFSVK": {"vendor": "procal", "irt": 76.57},
    "GASDFLSFAVK": {"vendor": "procal", "irt": 77.42},
    "GEDLDFIYVVK": {"vendor": "procal", "irt": 79.62},
    "VSSIFFDTFDK": {"vendor": "procal", "irt": 82.28},
    "SILDYVSLVEKK": {"vendor": "procal", "irt": 83.05},
    "VYGYELTSLFK": {"vendor": "procal", "irt": 87.89},
    "GGFFSFGDLTK": {"vendor": "procal", "irt": 88.04},
    "YDTAIDFGLFK": {"vendor": "procal", "irt": 89.4},
    "IVLFELEGITK": {"vendor": "procal", "irt": 94.97},
    "GIEDYYIFFAK": {"vendor": "procal", "irt": 95.37},
    "SILDYVSLVEK": {"vendor": "procal", "irt": 96.26},
    "AFSDEFSYFFK": {"vendor": "procal", "irt": 99.13},
    "AFLYEIIDIGK": {"vendor": "procal", "irt": 99.61},
}

SplitSet = Literal["Train", "Test", "Val"]

# Generated using {x:hash(x) for x in CONFIG.encoding_aa_order} once
# and then hard-coded
HASHDICT = {
    "A": 8990350376580739186,
    "C": -5648131828304525110,
    "D": 6043088297348140225,
    "E": 2424930106316864185,
    "F": 7046537624574876942,
    "G": 3340710540999258202,
    "H": 6743161139278114243,
    "I": -3034276714411840744,
    "K": -6360745720327592128,
    "L": -5980349674681488316,
    "M": -5782039407703521972,
    "N": -5469935875943994788,
    "P": -9131389159066742055,
    "Q": -3988780601193558504,
    "R": -961126793936120965,
    "S": 8601576106333056321,
    "T": -826347925826021181,
    "V": 6418718798924587169,
    "W": -3331112299842267173,
    "X": -7457703884378074688,
    "Y": 2606728663468607544,
}


def select_split(pep: str) -> SplitSet:
    """Assigns a peptide to a split set based on its sequence.

    It selects all iRT peptides to the 'Val' set.
    The rest of the peptides are hashed based on their stripped sequence (no mods).
    It is done on a semi-random basis

    Parameters
    ----------
    pep : str
        Peptide to assign to a split set

    Returns
    -------
    SplitSet: Split set to assign the peptide to.
        This is either one of "Train", "Test" or "Val"

    Exaples
    --------
    >>> select_split("AAA")
    'Train'
    >>> select_split("AAAK")
    'Test'
    >>> select_split("AAAKK")
    'Train'
    >>> select_split("AAAMTKK")
    'Train'

    """
    num_hash = sum(HASHDICT[x] for x in pep)
    num_hash = num_hash / 1e4
    num_hash = num_hash % 1
    assert 0 <= num_hash <= 1
    return _select_split(pep, num_hash)


# This is a hard-coded set of peptides that I use as landmarks to align
# My retention times but actively exclude in training.
def _select_split(pep: str, num_hash: float) -> SplitSet:
    """Assigns a peptide to a split set based on its sequence.

    It selects all iRT peptides to the 'Val' set.
    The rest of the peptides are hashed based on their stripped sequence (no mods).
    It is done on a semi-random basis, but is always the same between sessions.

    Parameters
    ----------
    pep : str
        Peptide to assign to a split set
    num_hash : float
        Hash of the peptide sequence

    Returns
    -------
        SplitSet: Split set to assign the peptide to.
            This is either one of "Train", "Test" or "Val"
    """
    in_landmark = pep in IRT_PEPTIDES
    if num_hash > 0.8 or in_landmark:
        return "Val"
    elif num_hash > 0.6:
        return "Test"
    else:
        return "Train"
