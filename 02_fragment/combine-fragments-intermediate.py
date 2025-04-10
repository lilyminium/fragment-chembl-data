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
from rdkit.Chem import rdChemReactions

def combine_fragment(fragment1, fragment2_smiles, reaction) -> str:
    """
    Combine two fragments with a single bond.
    """
    fragment2 = Chem.MolFromSmiles(fragment2_smiles)
    reactants = (fragment1, fragment2)
    
    products = [
        group[0]
        for group in reaction.RunReactants(reactants)
    ]
    return set([Chem.MolToSmiles(product) for product in products])


def combine_fragments_batch(
    index: int,
    fragment_smiles: list[str] = None,
    reaction = None
) -> set[str]:
    
    """
    Combine a single fragment with a list of other fragments.
    """
    combined = set()
    fragment1 = Chem.MolFromSmiles(fragment_smiles[index])
    for fragment2_smiles in fragment_smiles[index:]:
        combined|= combine_fragment(fragment1, fragment2_smiles, reaction)

    with open(f"intermediate/fragments-{index:08d}", "w") as f:
        f.write("\n".join(combined))


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
    "--n_processes",
    "-np",
    type=int,
    default=4,
    help="The number of processes to use for multiprocessing.",
    show_default=True,
)
def main(
    input_paths,
    intermediate_directory: str = "intermediate",
    n_processes: int = 4,
    n_heavy_atoms: int = 15,
):
    
    intermediate_directory = pathlib.Path(intermediate_directory)
    intermediate_directory.mkdir(parents=True, exist_ok=True)

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

    small_fragment_smiles = sorted(small_fragment_smiles, key=len, reverse=True)

    print(f"Filtered to {len(small_fragment_smiles)} smiles with {n_heavy_atoms} heavy atoms")

    reaction = rdChemReactions.ReactionFromSmarts("[*:1]-[#0].[*:3]-[#0]>>[*:1]-[*:3]")
    combiner = functools.partial(
        combine_fragments_batch,
        fragment_smiles=small_fragment_smiles,
        reaction=reaction
    )
    with multiprocessing.Pool(processes=n_processes) as pool:
        for fragment_set in tqdm.tqdm(
                pool.imap(combiner, range(len(small_fragment_smiles))),
                desc="Combining",
                total=len(small_fragment_smiles)
                ):
            pass



if __name__ == "__main__":
    main()
