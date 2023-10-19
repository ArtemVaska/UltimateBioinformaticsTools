from typing import List, Tuple, Dict, Union

from src import fastq_tools
from src import dna_rna_tools
from src import protein_tools


def filter_fastq_seqs(
        input_path: str,
        gc_bounds: Union[float, Tuple[float]] = (0, 100),
        length_bounds: Union[int, Tuple[int]] = (0, 2 ** 32),
        quality_threshold: int = 0,
        output_filename: str = None
) -> Dict[str, Tuple[str]]:
    """
    Filters provided sequences from .fastq file with specified parameters.

    Args:

    - input_path (str): name of file to filter without extension

    - gc_bounds (Union[float, Tuple[float]]): GC-content boundaries (from 0 to 100) within which filtered sequences must be included

    - length_bounds (Union[int, Tuple[int]]): length boundaries within which filtered sequences must be included

    - quality_threshold (int): reading quality value according to table phred+33, below which filtering will not be performed

    - output_filename (str): name of file where to save results without extension

    Return:

    - Dict[str, Tuple[str]]: a dictionary with filtered sequences
    """

    seqs = fastq_tools.read_file(input_path)

    if isinstance(gc_bounds, (int, float)):  # input check
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    if gc_bounds[1] > 100 or \
            (gc_bounds[0] > gc_bounds[1]) or \
            (length_bounds[0] > length_bounds[1]) or \
            not 0 <= quality_threshold <= 40:  # bounds check
        raise ValueError('The bounds are indicated incorrectly!')
    seqs_gc_and_len_in_bounds = (
        fastq_tools.is_gc_and_length_in_bounds(seqs, gc_bounds, length_bounds))
    quality_in_bounds = (
        fastq_tools.is_quality_in_bounds(seqs, quality_threshold))

    filtered_fastq_seqs = {}
    for key in seqs:
        if (seqs_gc_and_len_in_bounds[key][0] and
                seqs_gc_and_len_in_bounds[key][1] and
                quality_in_bounds[key]) is True:
            filtered_fastq_seqs[key] = seqs[key]

    if output_filename is None:
        output_filename = input_path
    fastq_tools.save_results(filtered_fastq_seqs, output_filename)

    return filtered_fastq_seqs


def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Accepts sequences and action to be done with them.

    Args:

    The last argument is the action performed on the sequences that are the other arguments

    Possible actions:

    - transcribe

    - complement

    - reverse

    - reverse_complement

    Return:

    - Union[str, List[str]]: sequence / sequences after action
    """

    *seqs, function = args
    functions = {
        'transcribe': dna_rna_tools.transcribe,
        'reverse': dna_rna_tools.reverse,
        'complement': dna_rna_tools.complement,
        'reverse_complement': dna_rna_tools.reverse_complement
    }
    result_list = []
    show_warning_message = False
    for seq in seqs:
        if ((set(seq) <= set('AUTGCautgc')) is False or
                'T' in seq.upper() and 'U' in seq.upper()):
            result_list.append(None)
            show_warning_message = True
        else:
            res = functions[function](seq)
            result_list.append(res)

    if show_warning_message:
        print('Some of your sequences have unreadable data!')

    if len(result_list) == 1:
        return result_list[0]
    return result_list


def run_ultimate_protein_tools(*args: Union[str, List[str]], **kwargs) -> list:
    """
    Accepts command and runs it on input data with parameters.

    Args:

    - args (str): the first argument is a command to do with sequences, which are other arguments

    - kwargs (str): various arguments for specific command

    Possible commands:

    - is_protein_valid

    - get_protein_rnas_number

    - get_length_of_protein

    - count_aa

    - get_fracture_of_aa

    Return:
    - list: a list with processed sequences
    """

    function, *seqs = args

    functions = {
        'is_protein_valid': protein_tools.is_protein_valid,
        'get_protein_rnas_number': protein_tools.get_protein_rnas_number,
        'get_length_of_protein': protein_tools.get_length_of_protein,
        'count_aa': protein_tools.count_aa,
        'get_fracture_of_aa': protein_tools.get_fracture_of_aa
    }
    result_list = []
    show_warning_message = False
    for seq in seqs:
        if not protein_tools.is_protein_valid(seq):
            result_list.append(None)
            show_warning_message = True
        else:
            res = functions[function](seq, **kwargs)
            result_list.append(res)
    if show_warning_message:
        print('Some of your sequences have mistakes!')
    if len(result_list) == 1:
        return result_list[0]
    return result_list
