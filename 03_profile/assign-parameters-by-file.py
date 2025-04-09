import multiprocessing
import pathlib
import logging
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

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

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
        entries.extend([
            {
                "smiles": smiles,
                "parameter_type": parameter_type,
                "parameter_id": parameter_id
            }
            for parameter_id in parameter_ids
        ])
    return entries


def label_single_file(
    input_file: str,
    output_directory: pathlib.Path,
    forcefield,
):
    output_file = pathlib.Path(input_file).stem + ".parquet"
    output_path =  output_directory / output_file

    if output_path.is_file():
        return

    with open(input_file, "r") as f:
        all_smiles = set([x.strip() for x in f.readlines()])

    logger.info(f"Loaded {len(all_smiles)} SMILES")

    entries = []
    for smiles in tqdm.tqdm(all_smiles):
        entries.extend(label_single_smiles(smiles, forcefield))
    
    # write table
    batch_table = pa.Table.from_pylist(entries)
    pq.write_table(batch_table, output_path)
    print(f"Wrote {len(entries)} entries to {output_path}")

def batch_label_files(
    all_files: list[str],
    forcefield_path: str,
    output_directory: str,
):

    forcefield = ForceField(forcefield_path)
    for ignore in ["vdW"]:
        forcefield.deregister_parameter_handler(ignore)
    output_directory = pathlib.Path(output_directory)
    for file_path in tqdm.tqdm(all_files):
        label_single_file(file_path, output_directory, forcefield)


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
    input_directory: str = "../02_fragment/intermediate",
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

    input_directory = pathlib.Path(input_directory)
    input_files = sorted(input_directory.glob("*"))[-20000:]


    print(f"Found {len(input_files)} files in {input_directory}")

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    with batch_distributed(
        input_files,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_label_files,
            forcefield_path=forcefield,
            output_directory=output_directory,
        ))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Updating entries",
        ):
            pass

    print("Done!")


if __name__ == "__main__":
    main()
