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


def read_gbk(input_gbk: str) -> Dict[str, List[str]]:
    """
    Reads /gene if possible or /locus_tag for each CDS in .gbk files.

    Additional function for select_genes_from_gbk_to_fasta.

    Args:

    - input_gbk (str): input filename without .gbk extension

    Return:

    - Dict[str, List[str]]: a dictionary with name of gene and its translation info
    """

    translations = {}
    with open(input_gbk + ".gbk") as gbk:
        translation = []
        cds_found = False
        translation_found = False
        gene_found = False
        # cds_counter = 0
        for line in gbk:
            if line.find("CDS") != -1:
                cds_found = True
                # cds_counter += 1
            elif cds_found:
                if not gene_found:
                    if line.find("/gene=") != -1:
                        gene = line[line.find("=") + 1 :].strip('"\n')
                    else:
                        gene = ""
                    gene_found = True

                if line.find("/locus_tag=") != -1:
                    locus_tag = line[line.find("=") + 1 :].strip('"\n')

                if line.find("/translation") != -1:
                    translation_found = True
                    name = gene if gene != "" else locus_tag
                    if line.count('"') == 2:
                        translations[name] = [
                            line[line.find("=") + 1 :].strip('"\n')
                        ]  # + [cds_counter]
                        translation = []
                        cds_found = False
                        translation_found = False
                        gene_found = False
                    else:
                        translation.extend(line[line.find('"') + 1 :].strip())
                elif translation_found:
                    if line.find('"') != -1:
                        translation.extend(line[21 : line.find('"')])
                        translations[name] = ["".join(translation)]  # + [cds_counter]
                        translation = []
                        cds_found = False
                        translation_found = False
                        gene_found = False
                    else:
                        translation.extend(line[21:].strip())
    return translations


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    # genes: List[str],
    # n_before: int = 1,
    # n_after: int = 1,
    output_fasta: str = None,
) -> Dict[str, List[str]]:
    """
    Writes found CDS' in .fasta extension.

    Args:

    - input_gbk (str): input filename without .gbk extension

    - output_fasta (str): output filename without .fasta extension

    Return:

    - Dict[str, List[str]]: a dictionary with name of gene and its translation info
    """

    translations = read_gbk(input_gbk)

    if output_fasta is None:
        warning_message = True
        output_fasta = input_gbk
    else:
        warning_message = False

    with open(output_fasta + ".fasta", mode="w") as output_fasta:
        for name, translation in translations.items():
            output_fasta.write(">" + name + "\n")
            output_fasta.write(translations[name][0] + "\n\n")

    if warning_message:
        print("Output_fasta wasn't provided.\n<input_gbk>.fasta created!")

    return translations
