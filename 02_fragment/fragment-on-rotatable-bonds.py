import logging
import multiprocessing
import pathlib
import typing
from collections import defaultdict

import click
from click_option_group import optgroup

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Recap
import tqdm


def canonical_smiles(rd_molecule: Chem.Mol) -> str:
    return Chem.MolToSmiles(
        Chem.AddHs(rd_molecule), isomericSmiles=False, canonical=True
    )


def fragment_single(parent_smiles: str):
    from openff.toolkit import Molecule

    # only keep longest part for things with ions
    parent_smiles = sorted(
        parent_smiles.split("."),
        key=len,
    )[-1]

    forbidden = {"Si", "BC", "Bc", "Al"}
    if any(el in parent_smiles for el in forbidden):
        return set()

    allowed_elements = {"Br", "C", "Cl", "F", "H", "I", "N", "O", "P", "S"}

    offmol = Molecule.from_smiles(parent_smiles, allow_undefined_stereo=True)
    rd_parent = offmol.to_rdkit()

    if any(
        rd_atom.GetNumRadicalElectrons() != 0
        or rd_atom.GetIsotope() != 0
        or rd_atom.GetSymbol() not in allowed_elements
        for rd_atom in rd_parent.GetAtoms()
    ):
        return set()
    
    bonds = offmol.find_rotatable_bonds()

    rdbonds = tuple([
        rd_parent.GetBondBetweenAtoms(bond.atom1_index, bond.atom2_index).GetIdx()
        for bond in bonds
    ])
    if len(rdbonds):
        fragmented = Chem.FragmentOnBonds(rd_parent, rdbonds)
    else:
        fragmented = rd_parent

    return set(Chem.MolToSmiles(fragmented).split("."))

@click.command()
@click.option(
    "--input",
    "-i",
    "input_paths",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    multiple=True,
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    default="fragments.smi",
    help="The file to the save the generated fragments to.",
    show_default=True,
)
@click.option(
    "--n_processes",
    "-np",
    type=int,
    default=4,
    help="The number of processes to use for multiprocessing.",
    show_default=True,
)
def main(
    input_paths,
    output_file,
    n_processes: int = 4,
):
    # Load in the molecules to fragment
    all_parent_smiles = set()

    for input_path in input_paths:
        with open(input_path, "r") as file:
            smiles = [x.strip() for x in file.readlines()]
            smiles = [
                sorted(x.split("."), key=len)[-1]
                for x in smiles
            ]
            all_parent_smiles |= set(smiles)

    all_parent_smiles = sorted(all_parent_smiles, key=len, reverse=True)
    print(f"Found {len(all_parent_smiles)} parent molecules")

    all_fragment_smiles = set()

    # Fragment the molecules with multiprocessing and tqdm
    with multiprocessing.Pool(processes=n_processes) as pool:
        all_sets = list(
            tqdm.tqdm(
                pool.imap(fragment_single, all_parent_smiles),
                total=len(all_parent_smiles),
                desc="Fragmenting",
            )
        )
        for fragment_set in all_sets:
            all_fragment_smiles |= fragment_set



    # for smi in tqdm.tqdm(all_parent_smiles):
    #     all_fragment_smiles |= fragment_single(smi)

    all_fragment_smiles = sorted(all_fragment_smiles, key=len, reverse=True)
    print(f"Found {len(all_fragment_smiles)} fragments")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as file:
        file.write("\n".join(all_fragment_smiles))
    print(f"Fragments written to {output_file}")



if __name__ == "__main__":
    main()
