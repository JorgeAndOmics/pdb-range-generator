# PDB Range Generator

## Description

PDB Range Generator is a Python script that allows users to extract a range of residues from PDB files for a protein of interest. The script downloads the amino acid sequence and PDB files from the UniProt and PDB websites, respectively, and selects the chain from each PDB file that shares the highest sequence identity to the amino acid sequence. It then extracts the desired range of residues from each chain and generates new PDB files with only the desired range of residues. The new PDB files retain the original ligands, solvent molecules, and other important information.

## Usage

To use the script, run the following command:

`python pdb-ranges-generator.py uniprot_code start_residue end_residue output_dir` 

| Argument | Description |
| -------- | ----------- |
| `uniprot_code` | The UniProt code for the protein of interest. |
| `start_residue` | The starting residue of the desired range. |
| `end_residue` | The ending residue of the desired range. |
| `output_dir` | The directory where the newly generated PDB files will be saved. |


## Usage Examples

Here are some examples of how to use the script:


| Command | Description |
| ------- | ----------- |
| `python extract_range_from_pdb.py P68871 10 20 output/` | Extract the range of residues 10-20 from the PDB files for the protein with UniProt code P68871, and save the newly generated PDB files in the `output/` directory. |
| `python extract_range_from_pdb.py Q9H3D4 50 100 results/` | Extract the range of residues 50-100 from the PDB files for the protein with UniProt code Q9H3D4, and save the newly generated PDB files in the `results/` directory. |


## Contribution

Contributions to this script are welcome! If you find a bug or have an idea for a new feature, please create an issue or submit a pull request.

## MIT License

The PDB Range Generator script is licensed under the MIT License. See the `LICENSE` file for more information.

## References

-   UniProt: [https://www.uniprot.org/](https://www.uniprot.org/)
-   Protein Data Bank: [https://www.rcsb.org/](https://www.rcsb.org/)
