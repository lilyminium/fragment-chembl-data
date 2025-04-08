"""
Combine single fragments with a single bond.
"""
import pathlib
import typing

import click
from click_option_group import optgroup

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


def batch_combine_fragments(
    indices: list[int],
    fragment_smiles: list[str] = None,
    reaction = None,
):
    """
    Combine a single fragment with a list of other fragments.
    """
    combined = set()
    
    for index in indices:
        fragment1 = Chem.MolFromSmiles(fragment_smiles[index])
        other_fragments = fragment_smiles[index:]
        
        for fragment2 in other_fragments:
            combined |= combine_fragment(fragment1, fragment2, reaction)
    
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
    "output_directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    help="The directory to save the generated fragments to.",
    show_default=True,
)
@click.option(
    "--n-heavy-atoms",
    "-n",
    help="The number of heavy atoms in the fragments to generate.",
    type=int,
    default=15,
    show_default=True,
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. Use -1 to request "
    "one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a particular worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("Cluster configuration", help="Options to configure cluster workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The SLURM queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that SLURM workers should run using.",
    type=str,
)
def main(
    input_paths: list[str],
    output_directory: str,
    n_heavy_atoms: int = 15,
    worker_type: typing.Literal["slurm", "local"] = "local",
    queue: str = "free",
    conda_environment: str = "ib-dev",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
    
):
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

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

    indices = list(range(len(small_fragment_smiles)))
    reaction = rdChemReactions.ReactionFromSmarts("[*:1]-[#0].[*:3]-[#0]>>[*:1]-[*:3]")

    with batch_distributed(
        indices,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_combine_fragments,
            fragment_smiles=small_fragment_smiles,
            reaction=reaction,
        ))
        for i, future in tqdm.tqdm(
            enumerate(
                distributed.as_completed(futures, raise_errors=False),
                0,
            ),
            total=len(futures),
            desc="Updating entries",
        ):
            batch = sorted(future.result(), key=len, reverse=True)
            output_file = output_directory / f"batch-{i:08d}.smi"
            with open(output_file, "w") as file:
                file.write("\n".join(batch))

            print(f"Wrote {len(batch)} to {output_file}")
    
    print("Done!")


if __name__ == "__main__":
    main()
