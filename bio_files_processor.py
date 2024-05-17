from dataclasses import dataclass


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> dict:
    """
    Creates oneline fasta from multiline fasta.
    :param input_fasta: input filename
    :param output_fasta: output filename; <input_filename>_oneline.fasta if not provided
    :return: a dictionary with id of sequence and its sequence
    """

    seqs = {}
    with open(input_fasta) as multiline_fasta:
        seq = []
        for line in multiline_fasta:
            line = line.strip()
            if line.startswith(">"):
                if len(seq) != 0:
                    seqs[seq_id] = "".join(seq)
                    seq = []
                seq_id = line
            else:
                seq.extend(line)
        seqs[seq_id] = "".join(seq)

    if output_fasta is None:
        output_fasta = input_fasta.split(".")[0] + "_oneline.fasta"

    with open(output_fasta, mode="w") as oneline_fasta:
        for seq_id, seq in seqs.items():
            oneline_fasta.write(f"{seq_id}\n")
            oneline_fasta.write(f"{seq}\n")

    return seqs


@dataclass
class FastaRecord:
    id_: str
    description: str
    seq: str

    def __repr__(self):
        header = f">{self.id_} {self.description}"
        return f"{header}\n{self.seq}"


class OpenFasta:
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
