from typing import Tuple, Union

import numpy as np
from Bio import SeqIO, SeqUtils


def filter_fastq(
    input_filename: str,
    output_filename: str = None,
    gc_bounds: Union[int | float, Tuple[int | float, int | float]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:

    if isinstance(gc_bounds, (int | float)):  # input check
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    if (
        gc_bounds[1] > 100
        or (gc_bounds[0] > gc_bounds[1])
        or (length_bounds[0] > length_bounds[1])
        or not 0 <= quality_threshold <= 40
    ):  # bounds check
        raise ValueError("The bounds are indicated incorrectly!")

    count_inp = 0
    for _ in SeqIO.parse(input_filename, "fastq"):
        count_inp += 1

    filtered_seqs = (
        record
        for record in SeqIO.parse(input_filename, "fastq")
        if np.mean(record.letter_annotations["phred_quality"]) >= quality_threshold
        and length_bounds[0] <= len(record.seq) <= length_bounds[1]
        and gc_bounds[0] <= SeqUtils.GC123(record)[0] <= gc_bounds[1]
    )

    if output_filename is None:
        output_filename = input_filename

    count_out = SeqIO.write(filtered_seqs, output_filename, "fastq")

    print(
        "Saved %i / %i sequences from %s in %s."
        % (count_out, count_inp, input_filename, output_filename)
    )


filter_fastq(
    "example.fastq",
    "result.fastq",
    gc_bounds=(50, 60),
    length_bounds=(220, 260),
    quality_threshold=10,
)


# def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
#     """
#     Accepts sequences and action to be done with them.
#
#     Args:
#
#     The last argument is the action performed on the sequences that are the other arguments
#
#     Possible actions:
#
#     - transcribe
#
#     - complement
#
#     - reverse
#
#     - reverse_complement
#
#     Return:
#
#     - Union[str, List[str]]: sequence / sequences after action
#     """
#
#     *seqs, function = args
#     functions = {
#         "transcribe": dna_rna_tools.transcribe,
#         "reverse": dna_rna_tools.reverse,
#         "complement": dna_rna_tools.complement,
#         "reverse_complement": dna_rna_tools.reverse_complement,
#     }
#     result_list = []
#     show_warning_message = False
#     for seq in seqs:
#         if (
#             (set(seq) <= set("AUTGCautgc")) is False
#             or "T" in seq.upper()
#             and "U" in seq.upper()
#         ):
#             result_list.append(None)
#             show_warning_message = True
#         else:
#             res = functions[function](seq)
#             result_list.append(res)
#
#     if show_warning_message:
#         print("Some of your sequences have unreadable data!")
#
#     if len(result_list) == 1:
#         return result_list[0]
#     return result_list
#
#
# def run_ultimate_protein_tools(*args: Union[str, List[str]], **kwargs) -> list:
#     """
#     Accepts command and runs it on input data with parameters.
#
#     Args:
#
#     - args (str): the first argument is a command to do with sequences, which are other arguments
#
#     - kwargs (str): various arguments for specific command
#
#     Possible commands:
#
#     - is_protein_valid
#
#     - get_protein_rnas_number
#
#     - get_length_of_protein
#
#     - count_aa
#
#     - get_fracture_of_aa
#
#     Return:
#     - list: a list with processed sequences
#     """
#
#     function, *seqs = args
#
#     functions = {
#         "is_protein_valid": protein_tools.is_protein_valid,
#         "get_protein_rnas_number": protein_tools.get_protein_rnas_number,
#         "get_length_of_protein": protein_tools.get_length_of_protein,
#         "count_aa": protein_tools.count_aa,
#         "get_fracture_of_aa": protein_tools.get_fracture_of_aa,
#     }
#     result_list = []
#     show_warning_message = False
#     for seq in seqs:
#         if not protein_tools.is_protein_valid(seq):
#             result_list.append(None)
#             show_warning_message = True
#         else:
#             res = functions[function](seq, **kwargs)
#             result_list.append(res)
#     if show_warning_message:
#         print("Some of your sequences have mistakes!")
#     if len(result_list) == 1:
#         return result_list[0]
#     return result_list
