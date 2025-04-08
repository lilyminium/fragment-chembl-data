import multiprocessing
import pathlib
import typing
import click
import tqdm

from click_option_group import optgroup

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pyarrow.compute as pc

from openff.toolkit import Molecule, ForceField

PARAMETER_TYPES = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"]

def label_single_smiles(smiles: str, forcefield):
    mol = Molecule.from_smiles(
        smiles,
        allow_undefined_stereo=True
    )
    labels = forcefield.label_molecules(mol.to_topology())[0]
    entries = []
    for parameter_type in PARAMETER_TYPES:
        values = labels[parameter_type]
        parameter_ids = sorted(set([
            parameter.id
            for parameter in values.values()
        ]))
        for parameter_id in parameter_ids:
            entry = {
                "smiles": smiles,
                "parameter_type": parameter_type,
                "parameter_id": parameter_id
            }
            entries.append(entry)
    return entries


def batch_label_smiles(
    all_smiles: list[str],
    forcefield_path: str,
):
    forcefield = ForceField(forcefield_path)
    entries = []
    for smiles in all_smiles:
        entries.extend(label_single_smiles(smiles, forcefield))
    return entries


@click.command()
@click.option(
    "--forcefield",
    "-f",
    "forcefield",
    default="openff_unconstrained-2.2.1.offxml",
    help="Force field to use for labeling",
)
@click.option(
    "--input",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="../02_fragment/combined",
    help="Directory containing input files",
)
@click.option(
    "--output",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="output",
    help="Directory to save output files",
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
    forcefield: str = "openff_unconstrained-2.2.1.offxml",
    input_file: str = "../02_fragment/output/fragments.smi",
    output_directory: str = "output",
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

    with open(input_file, "r") as f:
        all_smiles = set([x.strip() for x in f.readlines()])
    
    print(f"Found {len(all_smiles)} smiles")

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    start_index = 0

    try:
        print(f"Querying {output_directory}")
        existing = ds.dataset(output_directory)
        if existing.count_rows():
            seen_smiles = set(
                existing.to_table(columns=["smiles"]).to_pydict()["smiles"]
            )
            print(f"Found {len(seen_smiles)} existing smiles")
            all_smiles -= seen_smiles
            print(f"Filtered to {len(all_smiles)} new smiles")

            start_index = len(existing.files)

    except BaseException as e:
        print(e)

    all_smiles = sorted(all_smiles, key=len, reverse=True)[:200]

    with batch_distributed(
        all_smiles,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_label_smiles,
            forcefield_path=forcefield,
        ))
        for i, future in tqdm.tqdm(
            enumerate(
                distributed.as_completed(futures, raise_errors=False),
                start_index,
            ),
            total=len(futures),
            desc="Updating entries",
        ):
            batch = future.result()

            batch_table = pa.Table.from_pylist(batch)
            table_path = output_directory / f"batch-{i:06d}.parquet"

            pq.write_table(batch_table, table_path)
            print(f"Wrote {len(batch)} to {table_path}")

    print("Done!")


if __name__ == "__main__":
    main()
