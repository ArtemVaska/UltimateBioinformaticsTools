from dataclasses import dataclass


@dataclass
class FastaRecord:
    id_: str
    description: str
    seq: str

    def __repr__(self):
        header = f">{self.id_} {self.description}"
        return f"{header}\n{self.seq}"


class OpenFasta:
    """
    A context manager for reading FASTA files.

    Attributes:
    [x] path_to_fasta: the path to the FASTA file.
    [x] mode: the mode for opening the file, default is "r"

    Methods:
    [x] read_record(): reads the next FASTA record; analogue of readline()
    [x] read_records(): reads all FASTA records and returns a list of FastaRecord objects; analogue of readlines()
    """

    def __init__(self, path_to_fasta: str, mode: str = "r"):
        self.path_to_fasta = path_to_fasta
        self.mode = mode
        self.current_line = ""
        self.stop = False

    def __enter__(self):
        self.handler = open(self.path_to_fasta, self.mode)
        return self

    def __next__(self):
        if self.stop:
            raise StopIteration

        id_, description, seq = "", "", ""

        if self.current_line == "":
            self.current_line = next(self.handler)

        while True:
            line_ = self.current_line
            line_ = line_.strip()
            if line_.startswith(">"):
                if id_ == "":
                    line_ = line_.split()
                    id_, description = line_[0][1:], " ".join(line_[1:])
                else:
                    return FastaRecord(id_, description, seq)
            else:
                seq += line_

            try:
                self.current_line = next(self.handler)
            except StopIteration:
                self.stop = True
                return FastaRecord(id_, description, seq)

    def read_record(self):
        return self.__next__()

    def read_records(self):
        records = list()
        for record in self.__iter__():
            records.append(record)
        return records

    def __iter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.handler:
            self.handler.close()


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
    """
    Creates oneline fasta from multiline fasta.
    :param input_fasta: input filename
    :param output_fasta: output filename; <input_filename>_oneline.fasta if not provided
    """

    if output_fasta is None:
        output_fasta = input_fasta.split(".")[0] + "_oneline.fasta"

    with OpenFasta(input_fasta) as multiline_fasta:
        with open(output_fasta, "w") as oneline_fasta:
            for record in multiline_fasta:
                oneline_fasta.write(f"{record.__repr__()}\n")
