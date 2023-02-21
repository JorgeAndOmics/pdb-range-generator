import argparse
import os
import tempfile
import shutil
from Bio import SeqIO
from Bio.PDB import PDBList, PDBParser, Selection, Superimposer

def download_sequence(uniprot_code):
    """Downloads the amino acid sequence for a given UniProt code.

    Args:
        uniprot_code (str): The UniProt code for the protein of interest.

    Returns:
        Bio.Seq: The amino acid sequence for the given UniProt code.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_code}.fasta"
    response = urllib.request.urlopen(url)
    record = SeqIO.read(response, "fasta")
    return record.seq

def download_pdbs(uniprot_code, temp_dir):
    """Downloads the PDB files for a given UniProt code and converts them to .pdb files.

    Args:
        uniprot_code (str): The UniProt code for the protein of interest.
        temp_dir (str): The temporary directory where the PDB files will be stored.
    """
    pdb_list = PDBList()
    pdb_list.retrieve_pdb_file(uniprot_code, pdir=temp_dir, file_format="ent", overwrite=True)
    for ent_file in os.listdir(temp_dir):
        pdb_id = ent_file[3:7]
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(pdb_id, os.path.join(temp_dir, ent_file))
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(os.path.join(temp_dir, f"{pdb_id}.pdb"))

def pairwise_alignment(ref_seq, chain_seq):
    """Performs a pairwise local alignment between a reference sequence and a chain sequence.

    Args:
        ref_seq (Bio.Seq): The reference sequence (i.e., the amino acid sequence).
        chain_seq (str): The sequence of a PDB chain.

    Returns:
        tuple: A tuple containing the aligned amino acid sequence for the reference sequence and the chain sequence.
    """
    from Bio.pairwise2 import align
    alignments = align.localxs(ref_seq, chain_seq, -0.5, -0.1, gap_char="-", force_generic=False)
    best_alignment = max(alignments, key=lambda x: x[2])
    return best_alignment[0], best_alignment[1]

def select_best_chain(ref_seq, pdb_file):
    """Selects the chain from a PDB file that shares the highest sequence identity with a reference sequence.

    Args:
        ref_seq (Bio.Seq): The reference sequence (i.e., the amino acid sequence).
        pdb_file (str): The path to a PDB file.

    Returns:
        Bio.PDB.Chain.Chain: The best chain for the given PDB file.
    """
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_file, pdb_file)
    best_chain = None
    best_identity = 0
    for chain in structure.get_chains():
        chain_seq = "".join([residue.get_resname() for residue in chain.get_residues()])
        alignment_a, alignment_b = pairwise_alignment(ref_seq, chain_seq)
        identity = sum(a == b for a, b in zip(alignment_a, alignment_b)) / len(alignment_a)
        if identity > best_identity:
            best_chain = chain
            best_identity = identity
    return best_chain

def extract_range(pdb_file, start_residue, end_residue):
    """Extracts a range of residues from a PDB file.

    Args:
        pdb_file (str): The path to a PDB file.
        start_residue (int): The starting residue of the desired range.
        end_residue (int): The ending residue of the desired range.

    Returns:
        Bio.PDB.Structure.Structure: A new structure object containing only the desired range of residues.
    """
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_file, pdb_file)
    atom_list = Selection.unfold_entities(structure, "A")
    residues_to_remove = [residue for residue in structure.get_residues() if residue.id[1] < start_residue or residue.id[1] > end_residue]
    for residue in residues_to_remove:
        for atom in residue:
            atom_list.remove(atom)
    structure = Structure.Structure(pdb_file)
    structure.add(atom_list)
    return structure

def main(uniprot_code, start_residue, end_residue, output_dir):
    """Main function.

    Args:
        uniprot_code (str): The UniProt code for the protein of interest.
        start_residue (int): The starting residue of the desired range.
        end_residue (int): The ending residue of the desired range.
        output_dir (str): The directory where the newly generated PDB files will be saved.

    Returns:
        list: A list of the paths to the newly generated PDB files.
    """
    # Download the amino acid sequence
    sequence = download_sequence(uniprot_code)

    # Create a temporary directory to store PDB files
    temp_dir = tempfile.mkdtemp()

    # Download and convert PDB files to .pdb format
    download_pdbs(uniprot_code, temp_dir)

    # Select the best chain for each PDB file
    best_chains = []
    for pdb_file in os.listdir(temp_dir):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(temp_dir, pdb_file)
            best_chain = select_best_chain(sequence, pdb_path)
            best_chains.append(best_chain)

    # Extract the desired range from each best chain
    new_pdbs = []
    for i, best_chain in enumerate(best_chains):
        pdb_file = os.listdir(temp_dir)[i]
        pdb_id = pdb_file.split(".")[0]
        new_pdb = extract_range(best_chain, start_residue, end_residue)
        new_pdb_path = os.path.join(output_dir, f"{pdb_id}_{start_residue}-{end_residue}.pdb")
        PDBIO().set_structure(new_pdb)
        PDBIO().save(new_pdb_path)
        new_pdbs.append(new_pdb_path)

    # Clean up the temporary directory
    shutil.rmtree(temp_dir)

    return new_pdbs

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract a range of residues from PDB files.")
    parser.add_argument("uniprot_code", type=str, help="UniProt code")
    parser.add_argument("start_residue", type=int, help="Start residue")
    parser.add_argument("end_residue", type=int, help="End residue")
    parser.add_argument("output_dir", type=str, help="Output directory")
    args = parser.parse_args()

    # Run the main function
    new_pdbs = main(args.uniprot_code, args.start_residue, args.end_residue, args.output_dir)

    # Print the paths to the new PDB files
    for new_pdb in new_pdbs:
        print(new_pdb)
