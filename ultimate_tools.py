from typing import Tuple, Union
from abc import ABC, abstractmethod
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


# filter_fastq(
#     "example.fastq",
#     "result.fastq",
#     gc_bounds=(50, 60),
#     length_bounds=(220, 260),
#     quality_threshold=10,
# )


class BiologicalSequence(ABC):

    @abstractmethod
    def is_alphabet_correct(self):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, item):
        pass

    @abstractmethod
    def __repr__(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        self.seq = seq

    def is_alphabet_correct(self):
        return set(self.seq).issubset("ATGCU")

    def complement(self):  # TODO
        pass

    def gc_content(self):  # TODO
        pass

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, item):
        return self.seq[item]

    def __repr__(self):
        return self.seq


seq = NucleicAcidSequence("ATGC")

print(seq)
