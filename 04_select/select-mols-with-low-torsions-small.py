import click
import logging
import pathlib
import tqdm

import pyarrow as pa
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.compute as pc

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import SimDivFilters

from openff.toolkit import Molecule, ForceField

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    level=logging.INFO,
)

def parameter_tanimoto_distance(a: set[str], b: set[str]) -> float:
    """
    Calculate the Tanimoto distance between two sets of torsion parameters.
    """
    a = set(a)
    b = set(b)
    intersection = len(a.intersection(b))
    union = len(a.union(b))
    return 1 - intersection / union


def select_by_chemical_diversity(
    smiles: list[str],
    n_parameter_pool: int = 30,
) -> list[str]:
    """
    Select a subset of SMILES strings based on MorganFingerprint
    diversity using the MaxMinPicker.
    """
    smiles = sorted(smiles, key=len)
    rdmols = [Chem.MolFromSmiles(s) for s in smiles]
    fingerprints = [
        rdMolDescriptors.GetMorganFingerprintAsBitVect(rdmol, 2)
        for rdmol in rdmols
    ]
    mmp = SimDivFilters.MaxMinPicker()
    pick_size = min([n_parameter_pool, len(rdmols)])
    picked_indices = list(
        mmp.LazyBitVectorPick(
            fingerprints,
            poolSize=len(rdmols),
            pickSize=pick_size,
        )
    )
    return [smiles[i] for i in picked_indices]


def select_by_parameter_diversity(
    smiles: list[str],
    parameter_id: str,
    n_parameters: int = 10,
) -> list[str]:
    """
    Select a subset of SMILES strings based on Jaccard distance
    between torsions running through the same central bond.

    A li'l bit hacky.
    """
    ff = ForceField("openff_unconstrained-2.2.1.offxml")
    all_other_parameters = []
    unique_coupled_parameter_ids = set()
    ix_to_ignore = []
    for i, smi in enumerate(smiles):
        mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        indices = [
            i for i, value in labels.items()
            if value.id == parameter_id
        ]
        central_bonds = set()
        for ix in indices:
            # filter for being in ring
            i, j = ix[1:3]
            if not mol.get_bond_between(i, j).is_in_ring():
                central_bonds.add(tuple(sorted([i, j])))
        
        if not central_bonds:
            ix_to_ignore.append(i)
            
        other_parameter_ids = set([
            value.id for i, value in labels.items()
            if tuple(sorted(i[1:3])) in central_bonds
        ])
        all_other_parameters.append(other_parameter_ids)

        if central_bonds:
            unique_coupled_parameter_ids |= other_parameter_ids
    
    # remove those without central bonds
    smiles = [smi for i, smi in enumerate(smiles) if i not in ix_to_ignore]
    all_other_parameters = [
        all_other_parameters[i] for i in range(len(all_other_parameters))
        if i not in ix_to_ignore
    ]
   
    print(f"Unique coupled parameter ids: {unique_coupled_parameter_ids}")

    mmp = SimDivFilters.MaxMinPicker()
    dist_func = lambda i, j: parameter_tanimoto_distance(
        all_other_parameters[i],
        all_other_parameters[j],
    )
    pick_size = min([n_parameters, len(smiles)])
    picked_indices = list(
        mmp.LazyPick(
            dist_func,
            len(smiles), # pool size
            pick_size, # pick size
        )
    )
    for i in picked_indices:
        print(smiles[i], all_other_parameters[i])
    return [smiles[i] for i in picked_indices]


@click.command()
@click.option(
    "--input-directory",
    "-i",
    default="../03_profile/ff-parameters",
    help="Directory containing the Pyarrow dataset.",
)
@click.option(
    "--csv-file",
    "-c",
    default="all-matching-low-torsions.csv",
    help="Output CSV file.",
)
@click.option(
    "--smiles-file",
    "-s",
    default="selected-torsion-molecules.smi",
    help="Output SMILES file.",
)
@click.option(
    "--torsion-id",
    "-t",
    default="t126",
    help="Torsion ID to select.",
)
def main(
    input_directory="../03_profile/ff-parameters",
    csv_file="all-matching-low-torsions.csv",
    smiles_file="selected-torsion-molecules.smi",
    torsion_id: str = "t126",
    n_pool: int = 5000,
    n_parameter_pool: int = 250,
    n_output_parameters: int = 10,
):
    
    csv_file = pathlib.Path(csv_file)
    if not csv_file.is_file():

        dataset = ds.dataset(input_directory)
        expression = pc.field("parameter_id") == torsion_id
        torsion_subset = dataset.filter(expression)
        smiles = set(
            torsion_subset.to_table(columns=["smiles"]).to_pydict()["smiles"]
        )
        if "" in smiles:
            smiles.remove("")
        print(f"Found {len(smiles)} SMILES")

        # initial sort for length and take first 50000
        smiles = sorted(smiles, key=len)[:50000]
    
        expression2 = pc.field("smiles").isin(smiles)
        subset2 = torsion_subset.filter(expression2)
        df = subset2.to_table(columns=["parameter_id", "smiles"]).to_pandas()
        logger.info(f"Loaded {len(df)} rows from dataset")

        mws = []

        # this is also a filter for validation
        for smi in tqdm.tqdm(df.smiles.values):
            try:
                mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
            except:
                mw = 1e6
            else:
                mw = sum([atom.mass for atom in mol.atoms]).m
            mws.append(mw)

        df["mw"] = mws

        df.to_csv(csv_file)
        logger.info(f"Raw dataset saved to {csv_file}")

    df = pd.read_csv(csv_file, index_col=0)
    
    # sort by mw and take first 100
    df = df.sort_values("mw").head(n_pool)
    logger.info(f"Filtered dataset to {len(df)} rows")

    # select by chemical diversity
    parameter_smiles = select_by_chemical_diversity(
        df["smiles"].values,
        n_parameter_pool=n_parameter_pool,
    )
    logger.info(f"Selected {len(parameter_smiles)} SMILES by chemical diversity")

    # now sort by parameter diversity
    parameter_smiles = select_by_parameter_diversity(
        parameter_smiles,
        parameter_id=torsion_id,
        n_parameters=n_output_parameters,
    )
    logger.info(f"Selected {len(parameter_smiles)} SMILES by parameter diversity")
    output_smiles = [
        f"{smi} {torsion_id}"
        for smi in parameter_smiles
    ]

    with open(smiles_file, "w") as f:
        f.write("\n".join(output_smiles))
    logger.info(f"Selected SMILES saved to {smiles_file}")
    logger.info(f"Total number of selected SMILES: {len(output_smiles)}")


if __name__ == "__main__":
    main()
    

    
