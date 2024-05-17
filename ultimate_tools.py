import datetime
import os
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple, Union, Callable

import numpy as np
import requests
from Bio import SeqIO, SeqUtils
from bs4 import BeautifulSoup
from dotenv import load_dotenv

from genscan_scripts import parse_soup, URL, ORGANISMS, EXON_CUTOFFS


class WrongAlphabetError(ValueError):
    pass


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

        result = set(self.seq).issubset(self.alphabet)

        if not result:
            raise WrongAlphabetError("Alphabet of input sequence is not correct")

        return result

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, item) -> str:
        return self.seq[item]

    def __repr__(self) -> str:
        return self.seq


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq: str):
        super().__init__(seq)

    def gc_content(self) -> float:
        return round((self.seq.count("G") + self.seq.count("C")) / len(self.seq), 2)


class DNASequence(NucleicAcidSequence):
    alphabet = set("ATGC")

    def __init__(self, seq: str):
        super().__init__(seq)

    def complement(self):
        complemented_seq = self.seq.translate(str.maketrans("ATGCU", "TACGA"))
        return type(self)(complemented_seq)

    def transcribe(self):
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


def telegram_logger(chat_id: int) -> Callable:
    def decorator(func):
        load_dotenv("secure.env")
        tg_url = "https://api.telegram.org/bot"

        def wrapper(*args, **kwargs):
            start = time.time()
            try:
                func(*args, **kwargs)
            except Exception as error:
                message = f"Function `{func.__name__}` failed with an exception:\n\n`{type(error).__name__}: {error}`"
            else:
                end = time.time()
                execution_time = end - start

                time_obj = datetime.datetime.fromtimestamp(execution_time, datetime.UTC)
                if execution_time // 86400 < 1:
                    formatted_time = time_obj.strftime("%H:%M:%S.%f")
                else:
                    formatted_time = time_obj.strftime("%-j days, %H:%M:%S")

                message = f"Function `{func.__name__}` finished in `{formatted_time}`"

            requests.post(
                tg_url + os.getenv("TG_API_TOKEN") + "/sendMessage",
                params={
                    "chat_id": chat_id,
                    "parse_mode": "MarkdownV2",
                    "text": message,
                },
            )

        return wrapper

    return decorator


@dataclass
class GenscanOutput:
    status: int
    cds_list: str
    intron_list: str | list
    exon_list: str | list

    def __repr__(self):
        return (
            f"Status: {self.status}\n\n"
            f"Predicted peptide:\n{self.cds_list}\n\n"
            f"Predicted introns:\n{self.intron_list}\n\n"
            f"Predicted exons:\n{self.exon_list}"
        )


def run_genscan(
    sequence: str = None,
    sequence_file: str = None,
    organism: str = "Vertebrate",
    exon_cutoff: float = 0.01,
    sequence_name: str = "",
) -> GenscanOutput:
    if sequence is None and sequence_file is None:
        raise ValueError("Either sequence or sequence_file must be specified")

    if organism not in ORGANISMS:
        raise ValueError("Organism must be one of {}".format(ORGANISMS))

    if exon_cutoff not in EXON_CUTOFFS:
        raise ValueError("Exon_cutoff must be one of {}".format(EXON_CUTOFFS))

    data = {
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-p": "Predicted peptides only",
    }

    files = None
    if sequence_file is None:
        data["-s"] = sequence
    else:
        files = {"-u": open(sequence_file, "rb")}

    response = requests.post(URL, data=data, files=files)
    if response.status_code != 200:
        raise ValueError("Response status code is {}".format(response.status_code))

    soup = BeautifulSoup(response.text, "html.parser")
    predicted_peptide, exons, introns = parse_soup(soup)

    return GenscanOutput(response.status_code, predicted_peptide, introns, exons)
