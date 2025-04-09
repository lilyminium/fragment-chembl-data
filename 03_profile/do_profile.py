import pathlib
import time
import click
import tqdm

from yammbs.checkmol import analyze_functional_groups, ChemicalEnvironment
from openff.toolkit.topology import Molecule

import pyarrow as pa
import pyarrow.dataset as ds

import seaborn as sns
import matplotlib.pyplot as plt

@click.command()
@click.option(
    "--input", "-i",
    "input_smiles",
    default="chembl_35.smi", help="Input SMILES file")
@click.option(
    "--output", "-o",
    "output_path",
    default="output/chembl-35", help="Output path"
)
@click.option(
    "--images", "-p",
    "image_path",
    default="images-35", help="Path to save images"
)
def main(
    input_smiles: str,
    output_path: str,
    image_path: str = "images"
):
    now = time.time()
    print(f"Starting at {time.ctime(now)}")

    mols = []
    # load one-by-one for progress bar and throw out any failures or elements we don't want
    with open(input_smiles, "r") as f:
        contents = f.readlines()

    ALLOWED_ELEMENTS = ["C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"]

    for line in tqdm.tqdm(contents, desc="Loading molecules"):
        smi, chembl_id = line.strip().split()
        try:
            mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
            mol.name = chembl_id
            # filter out molecules with elements we don't want
            if any([a.symbol not in ALLOWED_ELEMENTS for a in mol.atoms]):
                continue
            mols.append(mol)
        except Exception as e:
            print(f"Failed to load {smi} with error {e}")
    print(f"Loaded {len(mols)} molecules")

    entries = []
    environments = [group.value for group in ChemicalEnvironment]
    base_environments = dict.fromkeys(
        ALLOWED_ELEMENTS + environments,
        False
    )

    for mol in tqdm.tqdm(mols):
        n_heavy_atoms = sum([1 for a in mol.atoms if a.atomic_number != 1])
        mw = sum([a.mass for a in mol.atoms]).m
        charge = mol.total_charge.m
        smi = mol.to_smiles()
        groups = [
            x.value for x in analyze_functional_groups(mol)
        ]
        entry = {
            "smiles": smi,
            "chembl_id": mol.name,
            "n_atoms": len(mol.atoms),
            "n_heavy_atoms": n_heavy_atoms,
            "mw": mw,
            "charge": charge,
            **base_environments
        }
        for atom in mol.atoms:
            entry[atom.symbol] = True
        for group in groups:
            entry[group] = True
        entries.append(entry)

    table = pa.Table.from_pylist(entries)
    ds.write_dataset(table, output_path)
    print(f"Wrote {len(entries)} entries to {output_path}")

    image_path = pathlib.Path(image_path)
    image_path.mkdir(parents=True, exist_ok=True)

    df = table.to_pandas()
    for cols in ["n_atoms", "n_heavy_atoms", "mw", "charge"]:
        print(f"Plotting {cols}")
        ax = sns.histplot(df[cols])
        ax.set_title(cols)
        plt.savefig(image_path / f"{cols}.png", dpi=300)
        print(f"Saved {cols}.png")
        plt.close()

    print(f"Finished at {time.ctime(time.time())}")
    print(f"Elapsed time: {time.time() - now} seconds")


if __name__ == "__main__":
    main()
