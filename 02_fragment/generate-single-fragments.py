"""
From fragments that may have multiple connection points, generate
fragments that have only one connection point.
"""


import multiprocessing
import pathlib

import click

from rdkit import Chem
import tqdm


def canonical_smiles(rd_molecule: Chem.Mol) -> str:
    return Chem.MolToSmiles(
        Chem.AddHs(rd_molecule), isomericSmiles=False, canonical=True
    )


def fragment_single(parent_smiles: str):
    rdmol = Chem.MolFromSmiles(parent_smiles)

    # get dummy atoms
    dummy_indices = [
        atom.GetIdx() for atom in rdmol.GetAtoms() if atom.GetAtomicNum() == "0"
    ]
    output_smiles = set()
    for idx in dummy_indices:
        copy = Chem.Mol(rdmol)
        for dummy_idx in dummy_indices:
            if dummy_idx == idx:
                continue
            copy.GetAtomWithIdx(dummy_idx).SetAtomicNum(1)
        output_smiles.add(Chem.MolToSmiles(copy))
    return output_smiles

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


    all_fragment_smiles = sorted(all_fragment_smiles, key=len, reverse=True)
    print(f"Found {len(all_fragment_smiles)} fragments")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as file:
        file.write("\n".join(all_fragment_smiles))
    print(f"Fragments written to {output_file}")



if __name__ == "__main__":
    main()
