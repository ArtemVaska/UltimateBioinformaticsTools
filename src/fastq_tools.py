from typing import List, Tuple, Dict

import os

ASCII_Q_SCORE = {
    '!': [33, 0],
    '"': [34, 1],
    '#': [35, 2],
    '$': [36, 3],
    '%': [37, 4],
    '&': [38, 5],
    "'": [39, 6],
    '(': [40, 7],
    ')': [41, 8],
    '*': [42, 9],
    '+': [43, 10],
    ',': [44, 11],
    '-': [45, 12],
    '.': [46, 13],
    '/': [47, 14],
    '0': [48, 15],
    '1': [49, 16],
    '2': [50, 17],
    '3': [51, 18],
    '4': [52, 19],
    '5': [53, 20],
    '6': [54, 21],
    '7': [55, 22],
    '8': [56, 23],
    '9': [57, 24],
    ':': [58, 25],
    ';': [59, 26],
    '<': [60, 27],
    '=': [61, 28],
    '>': [62, 29],
    '?': [63, 30],
    '@': [64, 31],
    'A': [65, 32],
    'B': [66, 33],
    'C': [67, 34],
    'D': [68, 35],
    'E': [69, 36],
    'F': [70, 37],
    'G': [71, 38],
    'H': [72, 39],
    'I': [73, 40]
    }

EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079804:1:SRR292678:1:1101:21885:21885':
        ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
         'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563':
        ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG',
         'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161':
        ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC',
         'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804:1:SRR292678:1:1101:47176:47176':
        ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG',
         'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079804:1:SRR292678:1:1101:149302:149302':
        ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA',
         '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
    '@SRX079804:1:SRR292678:1:1101:170868:170868':
        ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT',
         'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
    '@SRX079804:1:SRR292678:1:1101:171075:171075':
        ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT',
         'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
    '@SRX079804:1:SRR292678:1:1101:175500:175500':
        ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG',
         'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
    '@SRX079804:1:SRR292678:1:1101:190136:190136':
        ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT',
         'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
    '@SRX079804:1:SRR292678:1:1101:190845:190845':
        ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC',
         'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
    '@SRX079804:1:SRR292678:1:1101:198993:198993':
        ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA',
         '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
    '@SRX079804:1:SRR292678:1:1101:204480:204480':
        ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG',
         '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
    }


def calculate_gc_content(seq: str) -> float:
    """
    Calculates guanine-cytosine content in sequence.

    This is additional function for the filter_fastq_seqs function to work correctly

    Args:
    - seq (str): sequence in which to calculate GC-content

    Return:
    - float: GC-content in percentages
    """

    return (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100


def is_gc_and_length_in_bounds(
        seqs: Dict[str, Tuple[str]],
        gc_bounds: Tuple[float],
        length_bounds: Tuple[int]
) -> Dict[str, List[bool]]:
    """
    Checks whether sequences are within specified GC content and sequence length boundaries.

    This is additional function for the filter_fastq_seqs function to work correctly

    Args:
    - seqs (Dict[str, Tuple[str]]): a dictionary with FASTQ-file contents
    - gc_bounds (Tuple[float]): GC-content boundaries within which filtered sequences must be included
    - length_bounds (Tuple[int]): length boundaries within which filtered sequences must be included

    Return:
    - Dict[str, List[bool]]: a dictionary which contains the results of checking
        whether sequences in specified boundaries
    """

    seqs_gc_and_len_in_bounds = {}
    for key in seqs:
        seq = seqs[key][0]
        if gc_bounds[0] <= calculate_gc_content(seq) <= gc_bounds[1]:
            seqs_gc_and_len_in_bounds[key] = [True]
        else:
            seqs_gc_and_len_in_bounds[key] = [False]
        if length_bounds[0] <= len(seq) <= length_bounds[1]:
            seqs_gc_and_len_in_bounds[key].append(True)
        else:
            seqs_gc_and_len_in_bounds[key].append(False)
    return seqs_gc_and_len_in_bounds


def is_quality_in_bounds(
        seqs: Dict[str, Tuple[str]],
        quality_threshold: int = 0
) -> Dict[str, bool]:
    """
    Checks whether the average sequences read quality value is below the specified quality threshold.

    This is additional function for the filter_fastq_seqs function to work correctly

    Args:
    - seqs (Dict[str, Tuple[str]]): a dictionary with fastq-file contents
    - quality_threshold (int): sequences whose average read quality is lower than this will not be filtered

    Return:
    - Dict[str, bool]: a dictionary which contains the results of checking
        whether sequences are below the specified quality threshold
    """

    quality_in_bounds = {}
    for key in seqs:
        qual_seq = seqs[key][1]
        qual_list = []
        for value in qual_seq:
            qual_list.append(ASCII_Q_SCORE[value][1])
        mean_qual = sum(qual_list) / len(qual_list)
        if mean_qual >= quality_threshold:
            quality_in_bounds[key] = True
        else:
            quality_in_bounds[key] = False
    return quality_in_bounds


def read_file(
        input_path: str
) -> Dict[str, Tuple[str]]:
    """

    """

    seqs = {}
    with open(input_path) as fastq_file:
        for i, line in enumerate(fastq_file, start=1):
            line = line.strip()
            if i % 4 == 0:
                qual = line
                seqs[name] = (seq, qual)
            elif i % 2 == 0:
                seq = line
            elif i % 3 != 0 and line.startswith('@'):
                name = line
    return seqs


def save_results(
        filtered_fastq_seqs: Dict[str, Tuple[str]],
        output_filename: str
):
    """

    """

    dir_name = 'fastq_filtrator_results'
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    with open(os.path.join(dir_name, output_filename + '.fastq'), mode='w') as file:
        for name, (seq, qual) in filtered_fastq_seqs.items():
            file.write(name + '\n')
            file.write(seq + '\n')
            file.write('+' + '\n')
            file.write(qual + '\n')
    return
