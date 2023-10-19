from typing import Union, List

def convert_multiline_fasta_to_oneline(
        input_fasta: str,
        output_fasta: str = None
):
    """
    Creates

    Args:
        - input_fasta (str):
        - output_fasta (str):

    """

    seqs = {}
    with open(input_fasta + '.fasta') as multiline_fasta:
        seq = []
        for line in multiline_fasta:
            line = line.strip()
            if line.startswith('>'):
                if len(seq) != 0:
                    seqs[name] = ''.join(seq)
                    seq = []
                name = line
            else:
                seq.extend(line)
        seqs[name] = ''.join(seq)

    if output_fasta is None:
        warning_message = True
        output_fasta = input_fasta + '_oneline'
    else:
        warning_message = False

    with open(output_fasta + '.fasta', mode='w') as oneline_fasta:
        for name, seq in seqs.items():
            oneline_fasta.write(name + '\n')
            oneline_fasta.write(seq + '\n')

    if warning_message:
        print("Output_fasta wasn't provided.\n<input_fasta>_oneline created!")

    return seqs


def read_gbk(input_gbk: str):
    """

    """

    translations = {}
    with open(input_gbk + '.gbk') as gbk:
        translation = []
        cds_found = False
        translation_found = False
        gene_found = False
        for i, line in enumerate(gbk, start=1):
            if line.find('CDS') != -1:
                cds_found = True
            elif cds_found:
                if not gene_found:
                    if line.find('/gene=') != -1:
                        gene = line[line.find('=')+1:].strip('"\n')
                    else:
                        gene = ''
                    gene_found = True

                if line.find('/locus_tag=') != -1:
                    locus_tag = line[line.find('=')+1:].strip('"\n')

                if line.find('/translation') != -1:
                    translation_found = True
                    name = gene if gene != '' else locus_tag
                    if line.count('"') == 2:
                        translations[name] = line[line.find('=')+1:].strip('"')
                    else:
                        translation.extend(line[line.find('"')+1:])
                elif translation_found:
                    if line.find('"') != -1:
                        translation.extend(line[21:line.find('"')])
                        translations[name] = ''.join(translation)
                        translation = []
                        cds_found = False
                        translation_found = False
                        gene_found = False
                    else:
                        translation.extend(line[21:])
    return translations


def select_genes_from_gbk_to_fasta(
        input_gbk: str,
        # genes: Union[str, List[str]],
        # n_before: int,
        # n_after: int,
        output_fasta: str = None
):
    """

    """

    translations = read_gbk(input_gbk)

    # if output_fasta is None:
    #     warning_message = True
    #     output_fasta = input_gbk
    # else:
    #     warning_message = False
    #
    # with open(output_fasta + '.fasta', mode='w') as output_fasta:
    #     for name, translation in translations.items():
    #         output_fasta.write(name + '\n')
    #         output_fasta.write(translation + '\n')
    #
    # if warning_message:
    #     print("Output_fasta wasn't provided.\n<input_fasta>_oneline created!")

    return translations


translations = select_genes_from_gbk_to_fasta('example_gbk')

for key, value in translations.items():
    print(key, value, sep='===========')
