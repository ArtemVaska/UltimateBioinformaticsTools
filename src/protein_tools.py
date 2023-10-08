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


def get_protein_rnas_number(seq: str) -> int:
    """
    Gets number of all possible RNA's for a given protein.

    Args:
    - seq (str): sequence to count possible RNA's

    Return:
    - int: the number of possible RNA's for sequence
    """

    rnas_num = 1
    for amino_acid in seq:
        rnas_num *= len(RNA_AA_TABLE[amino_acid])
    return rnas_num


def get_length_of_protein(seq: str) -> int:
    """
    Calculates the length of a protein.

    Argument:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length
    """

    return len(seq)


def count_aa(
        seq: str,
        amino_acids: str = None
) -> Dict[str, int]:
    """
    Сounts the number of given or all amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence to count amino acids
    - aminoacids (str): which amino acids to count in sequence

    Return:
    - dict: a dictionary with amino acids and its count
    """

    aa_dict_count = {}
    if (amino_acids is None) or (amino_acids == ''):
        '''
        I added an additional condition for user-friendly experience.
        E.g., we can want to find specific aminoacid, look on result and then look on all aminoacids.
        Without this condition we have to delete keyword argument, but with it we can only make it empty.
        '''
        amino_acids = ''.join(set(seq))
    for aa in amino_acids:
        aa_dict_count[aa] = seq.count(aa)
    return aa_dict_count


def get_fracture_of_aa(
    seq: str,
    show_in_percentages: bool = False,
    amino_acids: str = None
) -> Dict[str, float]:
    """
    Calculates the fracture or percentage of amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Return:
    - dict: a dictionary with amino acids and its fracture or percentage with 2 or 4 digits after comma
    """

    if show_in_percentages:
        mult = 100
        round_var = 2
    else:
        mult = 1
        round_var = 4
    aa_dict_count = count_aa(seq, amino_acids=amino_acids)
    aa_dict_percent = {}
    len_of_protein = get_length_of_protein(seq)
    for aa, count in aa_dict_count.items():
        aa_dict_percent[aa] = round(count / len_of_protein * mult, round_var)
    return aa_dict_percent
