"""
Combine single fragments with a single bond.
"""
import itertools
import pathlib
import multiprocessing

import click
import tqdm

from rdkit import Chem

def combine_fragment(fragment_smiles: tuple[str, str]) -> str:
    """
    Combine two fragments with a single bond.
    """
    # Create an editable molecule from the two fragments
    fragment_1 = Chem.MolFromSmiles(fragment_smiles[0])
    fragment_2 = Chem.MolFromSmiles(fragment_smiles[1])
    rwmol = Chem.RWMol(
        Chem.CombineMols(fragment_1, fragment_2)
    )
    dummy_indices = [
        atom.GetIdx() for atom in rwmol.GetAtoms() if atom.GetAtomicNum() == 0
    ]
    if len(dummy_indices) != 2:
        raise ValueError("Expected exactly two dummy atoms to combine.")
    # Create a bond between the two dummy atoms
    rwmol.AddBond(dummy_indices[0], dummy_indices[1], Chem.BondType.SINGLE)
    
    # Sanitize the molecule
    Chem.SanitizeMol(rwmol)
    return Chem.MolToSmiles(rwmol)


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
    all_fragment_smiles = set()

    for input_path in input_paths:
        with open(input_path, "r") as file:
            smiles = [x.strip() for x in file.readlines()]
            all_fragment_smiles |= set(smiles)

    all_fragment_smiles = sorted(all_fragment_smiles, key=len, reverse=True)
    print(f"Found {len(all_fragment_smiles)} fragment molecules")

    all_combined_smiles = set()
    combinations = itertools.combinations_with_replacement(
        all_fragment_smiles, 2
    )

    # Fragment the molecules with multiprocessing and tqdm
    with multiprocessing.Pool(processes=n_processes) as pool:
        all_sets = list(
            tqdm.tqdm(
                pool.imap(combine_fragment, combinations),
                desc="Combining",
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
