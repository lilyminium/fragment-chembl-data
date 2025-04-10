# Fragmenting and combining ChEMBL molecules


Molecules were fragmented and combined in the following process:

- `fragment-on-rotatable-bonds.py`: Molecules were fragmented on all rotational bonds, using the OpenFF Toolkit `Molecule.find_rotatable_bonds()`. The output file is `fragments/chembl-35-fragments.smi`
- `generate-single-fragments.py`: Fragments were expanded to only have a single connection point per fragment. For each fragment with multiple dis/connection points, dummy atoms were replaced with hydrogens. The output file is `fragments/chembl-35-fragments-single.smi`
- `combine-fragments-intermediate.py`: Fragments with up to 15 heavy atoms were combined to generate unique molecules. Due to the combinatorial explosion of computational expense, intermediate files were saved to disk in `intermediate/` with the results. These files were not uploaded due to size limitations. 222075 files were generated