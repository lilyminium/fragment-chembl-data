"""
Combine single fragments with a single bond.
"""
import functools
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
    #fragment_1 = Chem.MolFromSmiles(fragment_smiles[0])
    #fragment_2 = Chem.MolFromSmiles(fragment_smiles[1])
    #rwmol = Chem.RWMol(
    #    Chem.CombineMols(fragment_1, fragment_2)
    #)
    rwmol = Chem.RWMol(Chem.MolFromSmiles(".".join(fragment_smiles)))
    dummy_indices = [
        atom.GetIdx() for atom in rwmol.GetAtoms() if atom.GetAtomicNum() == 0
    ]
    if len(dummy_indices) != 2:
        raise ValueError("Expected exactly two dummy atoms to combine.")
    # Create a bond between the two dummy atoms
    rwmol.AddBond(dummy_indices[0], dummy_indices[1], Chem.BondType.SINGLE)
    
    # Sanitize the molecule
    #Chem.SanitizeMol(rwmol)
    return Chem.MolToSmiles(rwmol)

def combine_fragments_batch(
    index: int,
    fragment_smiles: list[str],
) -> set[str]:
    
    """
    Combine a single fragment with a list of other fragments.
    """
    combined = set()
    fragment1 = fragment_smiles[index]
    other_fragments = fragment_smiles[index:]
    for fragment2 in other_fragments:
        combined.add(combine_fragment((fragment1, fragment2)))
    return combined


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
    n_heavy_atoms: int = 20,
):
    # Load in the molecules to fragment
    all_fragment_smiles = set()

    for input_path in input_paths:
        with open(input_path, "r") as file:
            smiles = [x.strip() for x in file.readlines()]
            all_fragment_smiles |= set(smiles)

    all_fragment_smiles = sorted(all_fragment_smiles, key=len, reverse=True)
    print(f"Found {len(all_fragment_smiles)} fragment molecules")
    small_fragment_smiles = set()
    for smi in tqdm.tqdm(all_fragment_smiles):
        rdmol = Chem.MolFromSmiles(smi)
        if rdmol.GetNumHeavyAtoms() <= n_heavy_atoms:
            small_fragment_smiles.add(smi)

    small_fragment_smiles = sorted(small_fragment_smiles, key=len, reverse=True)[:10]

    print(f"Filtered to {len(small_fragment_smiles)} smiles with {n_heavy_atoms} heavy atoms")

    all_combined_smiles = set()
    combiner = functools.partial(
        combine_fragments_batch,
        fragment_smiles=small_fragment_smiles,
    )
    # combinations = itertools.combinations_with_replacement(
    #     small_fragment_smiles, 2
    # )

    # Fragment the molecules with multiprocessing and tqdm
    with multiprocessing.Pool(processes=n_processes) as pool:
        all_sets = list(
            tqdm.tqdm(
                pool.imap(combiner, range(len(small_fragment_smiles))),
                desc="Combining",
            )
        )
        for fragment_set in all_sets:
            all_combined_smiles |= fragment_set


    all_combined_smiles = sorted(all_combined_smiles, key=len, reverse=True)
    print(f"Found {len(all_combined_smiles)} fragments")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as file:
        file.write("\n".join(all_combined_smiles))
    print(f"Fragments written to {output_file}")



if __name__ == "__main__":
    main()
