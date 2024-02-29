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
    """
    Filters provided sequences from .fastq file with specified parameters.

    Args:
        input_filename: name of input file (or path)
        output_filename: name of file with filtered seqs. If not provided creates filtered_<input_filename>
        gc_bounds: GC-content boundaries (from 0 to 100) to filter
        length_bounds: length of sequences boundaries to filter
        quality_threshold: phred33 quality threshold to filter

    Returns: creates the file with filtered sequences.

    Example:
        filter_fastq(
        "example.fastq",
        "result.fastq",
        gc_bounds=(50, 60),
        length_bounds=(220, 260),
        quality_threshold=10
        )
    """

    if isinstance(gc_bounds, (int | float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    if (
        gc_bounds[1] > 100
        or (gc_bounds[0] > gc_bounds[1])
        or (length_bounds[0] > length_bounds[1])
        or not 0 <= quality_threshold <= 40
    ):
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
        output_filename = "filtered_" + input_filename

    count_out = SeqIO.write(filtered_seqs, output_filename, "fastq")

    print(
        "Saved %i / %i sequences from %s in %s."
        % (count_out, count_inp, input_filename, output_filename)
    )


class AbstractBiologicalSequence(ABC):

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


class BiologicalSequence(AbstractBiologicalSequence):
    alphabet = None

    def __init__(self, seq: str):
        self.seq = seq

    def is_alphabet_correct(self) -> bool:
        if self.alphabet is None:
            raise NotImplementedError()
        return set(self.seq).issubset(self.alphabet)

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, item) -> str:
        return self.seq[item]

    def __repr__(self) -> str:
        return self.seq


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq: str):
        super().__init__(seq)

    def complement(self):
        if self.is_alphabet_correct():
            complemented_seq = self.seq.translate(str.maketrans("ATGCU", "TACGA"))
            return type(self)(complemented_seq)

    def gc_content(self) -> float:
        return round((self.seq.count("G") + self.seq.count("C")) / len(self.seq), 2)


class DNASequence(NucleicAcidSequence):
    alphabet = set("ATGC")

    def __init__(self, seq: str):
        super().__init__(seq)

    def transcribe(self):
        if self.is_alphabet_correct():
            transcribed_seq = self.seq.translate(str.maketrans("ATGC", "UACG"))
            return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence):
    alphabet = set("AUGC")

    def __init__(self, seq):
        super().__init__(seq)


class AminoAcidSequence(BiologicalSequence):
    alphabet = set("FLSY*CWPHQRIMTNKVADEG")

    def __init__(self, seq: str):
        super().__init__(seq)

    def count_aa(self) -> dict:
        amino_acids_dict = {}
        for aa in self.seq:
            amino_acids_dict[aa] = self.seq.count(aa)
        return amino_acids_dict
