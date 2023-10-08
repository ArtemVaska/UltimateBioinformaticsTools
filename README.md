# Almost-Ultimate-Bioinformatics-Tools
This is the repo for the 5th HomeWork of the BI Python 2023 course

### Overview

This project contains a `ultimate_tools.py` program, which implements the `run_ultimate_protein_tools`, `run_dna_rna_tools` and `filter_fastq_seqs` functions. These functions accept various biological sequences and the action that needs to be performed with them. 

### Installation

```python
import ultimate_tools as ut
```

Make sure the path to the directory with `ultimate_tools.py` is added to the PATH so that Python can find it when importing.

### Usage

1. **`run_ultimate_protein_tools`**

Accepts command and runs it on input data with parameters. 

    Args:
    - args (str): the first argument is a command to do with sequences, which are other arguments
    - kwargs (str): various arguments for specific command

    Possible commands:
    - is_protein_valid
    - get_protein_rnas_number
    - get_length_of_protein
    - count_aa
    - get_fracture_of_aa

    Return:
    - list: a list with processed sequences

In the first place in parentheses, you need to indicate the command from the possible ones, then the analyzed protein sequences in the form of strings, after which the named arguments necessary to execute the command, which can be viewed in the docstrings of the corresponding functions. The output will be a list of all given sequences.

1.1 `is_protein_valid`

Checks if protein is valid.

    Args:
    - seq (str): sequence to check

    Return:
    - bool: the result of the check

1.2 `get_protein_rnas_number`

Gets number of all possible RNA's for a given protein.

    Args:
    - seq (str): sequence to count possible RNA's

    Return:
    - int: the number of possible RNA's for sequence

1.3 `get_length_of_protein`

Calculates the length of a protein.

    Argument:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length

1.4 `count_aa`

Сounts the number of given or all amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence to count amino acids
    - aminoacids (str): which amino acids to count in sequence

    Return:
    - dict: a dictionary with amino acids and its count

1.5 `get_fracture_of_aa`

Calculates the fracture or percentage of amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Return:
    - dict: a dictionary with amino acids and its fracture or percentage with 2 or 4 digits after comma

2. **`run_dna_rna_tools`**

Accepts sequences and action to be done with them.

    Args:
    The last argument is the action performed on the sequences that are the other arguments

    Possible actions:
    - transcribe
    - complement
    - reverse
    - reverse_complement

    Return:
    - Union[str, List[str]]: sequence / sequences after action

First, in parentheses, you need to indicate all the analyzed DNA or RNA sequences in a string data type, and in the last place, the action from the listed possible ones that must be performed on all transferred nucleic sequences.

2.1 `transcribe`

Transcribes the sequence from DNA to RNA and vice versa.

    Args:
    - seq (str): a sequence to transcribe

    Return:
    - seq (str): transcribed sequence

2.2 `complement`

Returns a complementary sequence of DNA or RNA.

    Args:
    - seq (str): a sequence to complement

    Return:
    - seq (str): complemented sequence

2.3 `reverse`

Reverses the sequence.

    Args:
    - seq (str): a sequence to reverse

    Return:
    - seq (str): reversed sequence

2.4 `reverse_complement`

Returns a reversed complementary sequence of DNA or RNA.

    Args:
    - seq (str): a sequence to be reverse complemented

    Return:
    - seq (str): reverse complemented sequence

3. **`filter_fastq_seqs`**

Filters provided sequences with specified parameters.

    Args:
    - seqs (Dict[str, Tuple[str]]): a dictionary with FASTQ-file contents
    - gc_bounds (Union[float, Tuple[float]]): GC-content boundaries (from 0 to 100)
        within which filtered sequences must be included
    - length_bounds (Union[int, Tuple[int]]): length boundaries
        within which filtered sequences must be included
    - quality_threshold (int): reading quality value according to table phred+33,
        below which filtering will not be performed

    Return:
    - Dict[str, Tuple[str]]: a dictionary with filtered sequences

### Input of data

Functions `run_ultimate_protein_tools` and `run_dna_rna_tools` checks the input for mistakes and says if there are some of that, returning None for a specific sequence. 

Function `filter_fastq_seqs` assumes that the input is pre-processed data from the FASTQ file.

During each run of the first 2 functions, the user is required to enter a **protein sequence / sequences** or **DNA / RNA sequence / sequences** that must be processed using the procedures listed above.

The program involves the analysis of protein sequences consisting of <u>**20 canonical amino acids**</u> and nucleic acids consisting of <u>**5 nucleotides (A, T, G, C, U)**</u>.

```python
run_ultimate_protein_tools('get_protein_rnas_number', ['AAAAAAAAA', 'HJKASDKHSJAD'])
'Some of your sequences have mistakes!'
[262144, None]
```

### Examples

You can use the `EXAMPLE_FASTQ` constant for your tests!

```python
is_quality_in_bounds_mean_values(seqs=EXAMPLE_FASTQ)  # no parameteres = all sequences
{'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
  'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
 '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG',
  'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
 '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC',
  'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
 '@SRX079804:1:SRR292678:1:1101:47176:47176': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG',
  'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
 '@SRX079804:1:SRR292678:1:1101:149302:149302': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA',
  '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
 '@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT',
  'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
 '@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT',
  'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
 '@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG',
  'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
 '@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT',
  'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
 '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC',
  'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
 '@SRX079804:1:SRR292678:1:1101:198993:198993': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA',
  '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
 '@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG',
  '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')}

filter_fastq_seqs(fastq_tools.EXAMPLE_FASTQ, gc_bounds=(20, 50), length_bounds=75, quality_threshold=32)
{'@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG',
  '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')}

filter_fastq_seqs(seqs=EXAMPLE_FASTQ, length_bounds=(89, 89), quality_threshold=35)  # if you want to indicate a specific sequence size or GC composition, then indicate it twice as an interval in parentheses
{'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
  'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD')}

run_dna_rna_tools(fastq_tools.EXAMPLE_FASTQ['@SRX079804:1:SRR292678:1:1101:204480:204480'][0], 'reverse_complement')
'CGGGTCATTTTAAAGTGTGGCTAACTAGTAATATTCAAAGATGTCTTACTAGGAATGTTCAGGGGTGTCTCACT'

run_dna_rna_tools('ATG', 'transcribe')
'AUG'

run_dna_rna_tools('ATG', 'reverse')
'GTA'

run_dna_rna_tools('AtG', 'complement')
'TaC'

run_dna_rna_tools('ATG', 'aT', 'reverse')
['GTA', 'Ta']

run_ultimate_protein_tools('get_protein_rnas_number', ['AAAAAAAAA', 'HJKASDKHSJAD'])
'Some of your sequences have mistakes!'
[262144, None]

run_ultimate_protein_tools('get_fracture_of_aa', 'MMMfasdLA')
[{'M': 1.0},
 {'M': 1.0},
 {'M': 1.0},
 {'f': 1.0},
 {'a': 1.0},
 {'s': 1.0},
 {'d': 1.0},
 {'L': 1.0},
 {'A': 1.0}]

run_ultimate_protein_tools('is_protein_valid', 'MMMfasdLA', 'XXXDDDD')
'Some of your sequences have mistakes!'
[True, None]

run_ultimate_protein_tools('count_aa', 'MMMfasdLA')
{'M': 3, 'a': 1, 's': 1, 'f': 1, 'L': 1, 'A': 1, 'd': 1}
```

### Contacts
Artem Vasilev (artem_vasilev_01@list.ru)
