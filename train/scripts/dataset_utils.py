import re

from loguru import logger
from tqdm.auto import tqdm


def count_mods(peptide_sequences):
    mod_count = {"UNMODIFIED": 0}
    modification_regex = re.compile(r"\[([^\]]+)\]")
    if not any(modification_regex.search(peptide) for peptide in peptide_sequences):
        modification_regex = re.compile(r"\(([^\)]+)\)")

    for mod_peptide in peptide_sequences:
        mods = modification_regex.findall(mod_peptide)
        if not mods:
            mod_count["UNMODIFIED"] += 1
        for mod in mods:
            if mod not in mod_count:
                mod_count[mod] = 0
            mod_count[mod] += 1

    logger.info(f"Modifications: {mod_count}")

    return mod_count


def test_count_mods():
    # Count the number of modifications
    # _GM[Oxidation (M)]PVTAR = 1
    # _GM[Oxidation (M)]PVTAR[Oxidation (M)] = 2

    counts = count_mods(
        ["_GM[Oxidation (M)]PVTAR", "_GM[Oxidation (M)]PVTAR[Oxidation (M)]"],
    )
    assert counts["Oxidation (M)"] == 3

    counts = count_mods(["_(ac)AAAAAAGAGPEMVR_", "_(ac)AAAEVADTQLM(ox)LGVGLIEK_"])
    assert counts["ac"] == 2
    assert counts["ox"] == 1


def find_sequence_with_pairs(modified_peptides, stripped_peptides: list[str]):
    # Provided two lists of peptides of the same length,
    # one with modifications and one without, return a list of tuples
    # where the first element is the index to the modified peptide
    # and the second element is the index to the stripped peptide.

    # NOTE Pretty sure this function as a whole is prone to optimization.

    if "_" in modified_peptides[0]:
        stripped_peptides = [x.replace("_", "") for x in stripped_peptides]
        modified_peptides = [x.replace("_", "") for x in modified_peptides]

    modified_set = set(modified_peptides)

    out = []
    for i, (mp, sp) in tqdm(
        enumerate(zip(modified_peptides, stripped_peptides)),
        desc="Finding pairs",
        total=len(modified_peptides),
    ):
        # This makes sure the peptide is modified (not eq to stripped)
        # and that an unmodified version exists (stripped int he modified set).
        if ("(" in mp or "[" in mp) and mp != sp and sp in modified_set:
            out.append((i, modified_peptides.index(sp)))
    logger.info(f"Found {len(out)} pairs of modified and unmodified peptides")
    return out


def test_find_sequence_with_pairs():
    modified_peptides = [
        "_GM[Oxidation (M)]PVTAR",
        "_GM[Oxidation (M)]PVTAR[Oxidation (M)]",
    ]
    stripped_peptides = ["_GMPVTAR", "_GMPVTAR"]

    pairs = find_sequence_with_pairs(modified_peptides, stripped_peptides)
    assert pairs == []

    modified_peptides = ["_(ac)AAAAAAGAGPEMVR_", "_(ac)AAAEVADTQLM(ox)LGVGLIEK_"]
    stripped_peptides = ["AAAAAAGAGPEMVR", "AAAEVADTQLMLGVGLIEK"]

    pairs = find_sequence_with_pairs(modified_peptides, stripped_peptides)
    assert pairs == []

    modified_peptides = ["_(ac)AAAEVADTQLMLGVGLIEK", "AAAEVADTQLMLGVGLIEK"]
    stripped_peptides = ["AAAEVADTQLMLGVGLIEK", "AAAEVADTQLMLGVGLIEK"]

    pairs = find_sequence_with_pairs(modified_peptides, stripped_peptides)
    assert pairs == [(0, 1)]


def find_mod_pairs(modified_peptides, stripped_peptides):
    logger.info("Finding pairs of modifications")
    pair_indices = find_sequence_with_pairs(modified_peptides, stripped_peptides)
    counts = count_mods(modified_peptides)

    out = {}
    for mod in counts:
        out[mod] = []
        for i, j in pair_indices:
            if mod in modified_peptides[i]:
                out[mod].append((i, j))

    counts_per_mod = {mod: len(pairs) for mod, pairs in out.items()}
    logger.info(f"Found {counts_per_mod} pairs of modifications")
    return out


def test_find_mod_pairs():
    modified_peptides = ["_(ac)AAAEVADTQLMLGVGLIEK", "AAAEVADTQLMLGVGLIEK"]
    stripped_peptides = ["AAAEVADTQLMLGVGLIEK", "AAAEVADTQLMLGVGLIEK"]

    pairs_dict = find_mod_pairs(modified_peptides, stripped_peptides)
    assert pairs_dict["ac"] == [(0, 1)]

    modified_peptides = ["_GMPVTAR", "_GM[Oxidation (M)]PVTAR[Oxidation (M)]"]
    stripped_peptides = ["_GMPVTAR", "_GMPVTAR"]

    pairs_dict = find_mod_pairs(modified_peptides, stripped_peptides)
    assert pairs_dict["Oxidation (M)"] == [(1, 0)]


def write_simple_toml(file, data: dict[str, float]):
    with open(file, "w+") as f:
        for key, value in data.items():
            f.write(f"{key} = {value}\n")


# TODO abstract out evaluation logic.
