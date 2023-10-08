ALPHABET_FOR_TRANSCRIBE = {
    'T': 'U', 't': 'u',
    'U': 'T', 'u': 't'
    }

ALPHABET_FOR_COMPLEMENT = {
    'A': 'T', 'a': 't',
    'G': 'C', 'g': 'c',
    'T': 'A', 't': 'a',
    'C': 'G', 'c': 'g',
    'U': 'A', 'u': 'a'
    }


def transcribe(seq: str) -> str:
    """
    Transcribes the sequence from DNA to RNA and vice versa.

    Args:
    - seq (str): a sequence to transcribe

    Return:
    - seq (str): transcribed sequence
    """

    return ''.join([ALPHABET_FOR_TRANSCRIBE[nucleotide] if nucleotide in ALPHABET_FOR_TRANSCRIBE else
                    nucleotide for nucleotide in seq])


def complement(seq: str) -> str:
    """
    Returns a complementary sequence of DNA or RNA.

    Args:
    - seq (str): a sequence to complement

    Return:
    - seq (str): complemented sequence
    """

    return ''.join([ALPHABET_FOR_COMPLEMENT[nucleotide] for nucleotide in seq])


def reverse(seq: str) -> str:
    """
    Reverses the sequence.

    Args:
    - seq (str): a sequence to reverse

    Return:
    - seq (str): reversed sequence
    """

    return seq[::-1]


def reverse_complement(seq: str) -> str:
    """
    Returns a reversed complementary sequence of DNA or RNA.

    Args:
    - seq (str): a sequence to be reverse complemented

    Return:
    - seq (str): reverse complemented sequence
    """

    return reverse(complement(seq))
