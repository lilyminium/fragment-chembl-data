import time
import click

import rdkit
from rdkit import Chem

@click.command()
@click.option("--input_file", "-i", default="input/chembl_35.sdf", help="Input file")
@click.option("--output_file", "-o", default="output/chembl_35.smi", help="Output file")
@click.option("--n_processes", "-n", default=1, help="Number of processes to use")
def main(
    input_file: str = "input/chembl_35.sdf",
    output_file: str = "output/chembl.smi",
    n_processes: int = 1
):
    now = time.time()
    print(f"Starting at {time.ctime(now)}")

    print(f"RDKIT version: {rdkit.__version__}")

    mol_supplier = Chem.MultithreadedSDMolSupplier(input_file, n_processes)
    all_lines = []
    for mol in mol_supplier:
        if mol is not None:
            props = mol.GetPropsAsDict()
            all_lines.append(
                f"{Chem.MolToSmiles(mol)} {props['chembl_id']}"
            )
    with open(output_file, "w") as f:
        f.write("\n".join(all_lines))
    print(f"Converted {len(all_lines)} molecules to SMILES")

    print(f"Finished at {time.ctime(time.time())}")
    print(f"Elapsed time: {time.time() - now} seconds")

if __name__ == "__main__":
    main()
