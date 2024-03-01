"""
Uses libraries and packages from BioPython
"""
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB import Superimposer
from Bio.Data.PDBData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio import pairwise2
from Bio.Align import substitution_matrices

backbone_atoms = ["N","CA","C","O"]

def get_str_obj(pdb1,pdb2):
    # get the first model and first chain
    chain_lst1 = []
    chain_lst2 = []

    for chain in pdb1[0]:
        chainid = chain.get_id()
        chain_lst1.append(chainid)

    for chain in pdb2[0]:
        chainid = chain.get_id()
        chain_lst2.append(chainid)

    # define ref and mobile strs
    reference = pdb1[0][chain_lst1[0]]
    mobile = pdb2[0][chain_lst2[0]]

    return reference,mobile

def res_range(start_id,end_id):
    """
    Select what residues numbers you wish to align
    and put them in a list
    """
    atoms_to_be_aligned = range(int(start_id), int(end_id) + 1)
    #print(atoms_to_be_aligned)
    return(atoms_to_be_aligned)

def align_sequences(structA, structB):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    def _get_pdb_sequence(structure):
        """
        Retrieves the AA sequence from a PDB structure.
        """

        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, "X"))
        seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
        return seq

    resseq_A = _get_pdb_sequence(structA)
    resseq_B = _get_pdb_sequence(structB)

    sequence_A = "".join([i[1] for i in resseq_A])
    sequence_B = "".join([i[1] for i in resseq_B])

    alns = pairwise2.align.globalds(
        sequence_A,
        sequence_B,
        substitution_matrices.load("BLOSUM62"),
        one_alignment_only=True,
        open=-10.0,
        extend=-0.5,
        penalize_end_gaps=(False, False),
    )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == "-":
            if aa_aln_B != "-":
                aa_i_B += 1
        elif aa_aln_B == "-":
            if aa_aln_A != "-":
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
            aa_i_A += 1
            aa_i_B += 1

    return mapping

def align(pdb1: Structure, pdb2: Structure):
    """
    Aligns a model structure onto a native structure
    based on sequence alignment

    for whole chains

    """
    # define the reference and mobile str objects for the calculation
    reference, mobile = get_str_obj(pdb1,pdb2)

    # Align sequences to get mapping between residues
    mapping = align_sequences(reference, mobile)

    # residue list for subsequent calculation
    refe_atom_list, mobi_atom_list = [], []
    
    for refe_res in mapping:
        for atom in backbone_atoms:
            if atom in reference[refe_res] and atom in mobile[mapping[refe_res]]:
                refe_atom_list.append(reference[refe_res][atom])
                mobi_atom_list.append(mobile[mapping[refe_res]][atom])
    
    # Superimpose matching residues
    si = Superimposer()
    si.set_atoms(refe_atom_list, mobi_atom_list)
    si.apply(mobile.get_atoms())
    
    return si

def align_partial(pdb1: Structure, pdb2: Structure, res_range1:list, res_range2:list):
    """
    Aligns a model structure onto a native structure
    based on sequence alignment

    for partial regions

    """
    # define the reference and mobile str objects for the calculation
    reference, mobile = get_str_obj (pdb1,pdb2)

    # Align sequences to get mapping between residues
    mapping = align_sequences(reference, mobile)

    # residue list for subsequent calculation
    refe_atom_list, mobi_atom_list = [], []
    
    for refe_res in mapping:
        if refe_res in res_range1 and mapping[refe_res] in res_range2:
            for atom in backbone_atoms:
                if atom in reference[refe_res] and atom in mobile[mapping[refe_res]]:
                    refe_atom_list.append(reference[refe_res][atom])
                    mobi_atom_list.append(mobile[mapping[refe_res]][atom])
    
    if len(mobi_atom_list)==0:
        return 20
    if len(refe_atom_list)==0:
        return 20

    # Superimpose matching residues
    si = Superimposer()
    si.set_atoms(refe_atom_list, mobi_atom_list)
    si.apply(mobile.get_atoms())
    
    return si


def cal_rmsd(file1="/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data/clean_chains/2oug_C.pdb", file2="/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data/clean_chains/6c6s_D.pdb"):
    # define the PDB parser
    p = PDBParser(QUIET=True)
    # get the structure objects
    pdb1 = p.get_structure("reference", file1)
    pdb2  = p.get_structure("mobile", file2)
    si = align(pdb1, pdb2)
    if si==20:
        return 20
    else:
        return(si.rms)


def cal_rmsd_partial(file1="//Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data/clean_chains/2oug_C.pdb", file2="/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data/clean_chains/6c6s_D.pdb", start1=115, end1=162, start2=112, end2=162):
    # define the PDB parser
    p = PDBParser(QUIET=True)
    # get the structure objects
    pdb1 = p.get_structure("reference", file1)
    pdb2  = p.get_structure("mobile", file2)
    res_range1 = res_range(start1,end1)
    res_range2 = res_range(start2,end2)
    si = align_partial(pdb1,pdb2,res_range1,res_range2)
    if si == 20:
        return 20
    else:
        return(si.rms)

if __name__ == "__main__": 
    rmsd = cal_rmsd()
    print(f"Whole RMSD: {rmsd:4.2f}")
    partial_rmsd = cal_rmsd_partial()
    print(f"Partial RMSD: {partial_rmsd:4.2f}")
