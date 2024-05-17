from typing import List, Dict


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
