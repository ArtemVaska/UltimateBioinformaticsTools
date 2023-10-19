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


convert_multiline_fasta_to_oneline('example_multiline_fasta')
