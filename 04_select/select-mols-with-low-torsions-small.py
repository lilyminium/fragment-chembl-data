import click
import logging
import tqdm

import pyarrow as pa
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
    picked_indices = list(
        mmp.LazyBitVectorPick(
            fingerprints,
            poolSize=len(rdmols),
            pickSize=n_parameter_pool,
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
    ff = ForceField("openff-2.2.1.offxml")
    all_other_parameters = []
    for smi in smiles:
        mol = Molecule.from_smiles(smi)
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        indices = [
            i for i, value in labels.items()
            if value.id == parameter_id
        ]
        central_bond = sorted(indices[1:3])
        other_parameter_ids = set([
            value.id for i, value in labels.items()
            if i[1:3] == central_bond
        ])
        all_other_parameters.append(other_parameter_ids)

    mmp = SimDivFilters.MaxMinPicker()
    dist_func = lambda i, j: parameter_tanimoto_distance(
        all_other_parameters[i],
        all_other_parameters[j],
    )
    picked_indices = list(
        mmp.LazyPick(
            dist_func,
            pool_size=len(smiles),
            pick_size=n_parameters,
        )
    )
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
    n_pool: int = 100,
    n_parameter_pool: int = 30,
    n_output_parameters: int = 10,
):
    dataset = ds.dataset(input_directory)
    expression = pc.field("parameter_id") == torsion_id
    torsion_subset = dataset.filter(expression)
    smiles = set(
        torsion_subset.to_table(columns=["smiles"]).to_pydict()["smiles"]
    )
    if "" in smiles:
        smiles.remove("")

    # initial sort for length and take first 5000
    smiles = sorted(smiles, key=len)[:5000]
    
    expression2 = pc.field("smiles").isin(smiles)
    subset2 = torsion_subset.filter(expression2)
    df = subset2.to_table(columns=["parameter_id", "smiles"]).to_pandas()
    logger.info(f"Loaded {len(df)} rows from dataset")

    # <= 10 torsions in central bonds across all QCArchive
    # TORSION_IDS = [
    #     't126', 't128', 't8', 't114', 't112', 't164', 't18b', 't81', 't141b',
    #     't138a', 't141a', 't113', 't167', 't137', 't7', 't87a', 't136', 't102',
    #     't33', 't141', 't89', 't103', 't31a', 't165', 't141c', 't55', 't54',
    #     't73', 't49', 't42a', 't101', 't88', 't158', 't12', 't154', 't129',
    #     't30'
    # ]
    # logger.info(f"Filtering dataset for {len(TORSION_IDS)} torsion IDs")

    # dfs = []
    # for torsion_id in tqdm.tqdm(TORSION_IDS):
    #     expression = pc.field("parameter_id") == torsion_id
    #     torsion_subset = dataset.filter(expression)
    #     smiles = torsion_subset.to_table(columns=["smiles"]).to_pydict()["smiles"]
    #     # initial sort for length and take first 5000
    #     smiles = sorted(smiles, key=len)[:5000]
        
    #     expression2 = pc.field("smiles").isin(smiles)
    #     subset2 = torsion_subset.filter(expression2)
    #     dfs.append(
    #         subset2.to_table(columns=["parameter_id", "smiles"]).to_pandas()
    #     )
    # df = pd.concat(dfs)

    #expression = pc.field("parameter_id").isin(TORSION_IDS)
    #subset = dataset.filter(expression)

    #df = subset.to_table().to_pandas()
    # logger.info(f"Loaded {len(df)} rows from dataset")

    # add mw
    mws = []

    # this is also a filter for validation
    for smi in tqdm.tqdm(df.smiles.values):
        try:
            mol = Molecule.from_smiles(smi)
        except:
            mw = 1e6
        else:
            mw = sum([atom.mass for atom in mol.atoms]).m
        mws.append(mw)

    df["mw"] = mw

    df.to_csv(csv_file)
    logger.info(f"Raw dataset saved to {csv_file}")

    
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
    

    
