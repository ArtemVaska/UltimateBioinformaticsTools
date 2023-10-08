from typing import Dict

AMINOACID_DICT = {
    'A': 'Alanine', 'a': 'alanine',
    'C': 'Cysteine', 'c': 'cysteine',
    'D': 'Aspartic acid', 'd': 'aspartic acid',
    'E': 'Glutamic acid', 'e': 'glutamic acid',
    'F': 'Phenylalanine', 'f': 'Phenylalanine',
    'G': 'Glycine', 'g': 'glycine',
    'H': 'Histidine', 'h': 'histidine',
    'I': 'Isoleucine', 'i': 'isoleucine',
    'K': 'Lysine', 'k': 'lysine',
    'L': 'Leucine', 'l': 'leucine',
    'M': 'Methionine', 'm': 'methionine',
    'N': 'Asparagine', 'n': 'asparagine',
    'P': 'Proline', 'p': 'proline',
    'Q': 'Glutamine', 'q': 'glutamine',
    'R': 'Arginine', 'r': 'arginine',
    'S': 'Serine',  's': 'serine',
    'T': 'Threonine', 't': 'threonine',
    'V': 'Valine', 'v': 'valine',
    'W': 'Tryptophan', 'w': 'tryptophan',
    'Y': 'Tyrosine', 'y': 'tyrosine'
    }

RNA_AA_TABLE = {
    'F': ['UUU', 'UUC'], 'f': ['uuu', 'uuc'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'l': ['uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 's': ['ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'],
    'Y': ['UAU', 'UAC'], 'y': ['uau', 'uac'],
    '*': ['UAA', 'UAG', 'UGA', 'uaa', 'uag', 'uga'],
    'C': ['UGU', 'UGC'], 'c': ['ugu', 'ugc'],
    'W': ['UGG'], 'w': ['ugg'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'p': ['ccu', 'ccc', 'cca', 'ccg'],
    'H': ['CAU', 'CAC'], 'h': ['cau', 'cac'],
    'Q': ['CAA', 'CAG'], 'q': ['caa', 'cag'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'r': ['cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
    'I': ['AUU', 'AUC', 'AUA'], 'i': ['auu', 'auc', 'aua'],
    'M': ['AUG'], 'm': ['aug'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'], 't': ['acu', 'acc', 'aca', 'acg'],
    'N': ['AAU', 'AAC'], 'n': ['aau', 'aac'],
    'K': ['AAA', 'AAG'], 'k': ['aaa', 'aag'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'v': ['guu', 'guc', 'gua', 'gug'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'a': ['gcu', 'gcc', 'gca', 'gcg'],
    'D': ['GAU', 'GAC'], 'd': ['gau', 'gac'],
    'E': ['GAA', 'GAG'], 'e': ['gaa', 'gag'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'], 'g': ['ggu', 'ggc', 'gga', 'ggg']
    }


def is_protein_valid(seq: str) -> bool:
    """
    Checks if protein is valid.

    Args:
    - seq (str): sequence to check

    Return:
    - bool: the result of the check
    """

    if set(seq).issubset(RNA_AA_TABLE):
        return True
    return False
