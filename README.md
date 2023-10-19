# Almost-Ultimate-Bioinformatics-Tools
This is the repo for the 6th HomeWork of the BI Python 2023 course

### Overview

This project contains a `ultimate_tools.py` program, which implements the `run_ultimate_protein_tools`, `run_dna_rna_tools` and `filter_fastq_seqs` functions. These functions accept various biological sequences and the action that needs to be performed with them. 

Also it contains `bio_file_processor.py` program, which provides next functions `convert_multiline_fasta_to_oneline`, `select_genes_from_gbk_to_fasta`. Unfortunately, second function for now only reads .gbk file and creates a new file in .fasta format with name of genes and its translation information.

### Installation

```python
import ultimate_tools as ut

import bio_files_processor as bfp
```

Make sure the path to the directory with `ultimate_tools.py` and `bio_files_processor` is added to the PATH so that Python can find it when importing.

All files you work with must be in the directory with the corresponding scripts.

### Usage

#### I. `ultimate_tools.py`

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

Filters provided sequences from .fastq file with specified parameters. Creates `fastq_filtrator_results` dir if it doesn't exist and save the result of filtration to this dir.

    Args:

    - input_path (str): name of file to filter without extension
    - gc_bounds (Union[float, Tuple[float]]): GC-content boundaries (from 0 to 100) within which filtered sequences must be included
    - length_bounds (Union[int, Tuple[int]]): length boundaries within which filtered sequences must be included
    - quality_threshold (int): reading quality value according to table phred+33, below which filtering will not be performed
    - output_filename (str): name of file where to save results without extension

    Return:

    - Dict[str, Tuple[str]]: a dictionary with filtered sequences

#### II. `bio_file_processor.py`

1. **`convert_multiline_fasta_to_oneline`**

Creates oneline fasta from multiline fasta.

    Args:

    - input_fasta (str): input filename without .fasta extension
    - output_fasta (str): output filename without .fasta extension

    Return:

    - Dict[str, str]: a dictionary with id of sequence and its sequence

2. **`read_gbk`**

Reads /gene if possible or /locus_tag for each CDS in .gbk files. Additional function for select_genes_from_gbk_to_fasta.

    Args:

    - input_path (str): name of file to filter without extension
    - gc_bounds (Union[float, Tuple[float]]): GC-content boundaries (from 0 to 100) within which filtered sequences must be included
    - length_bounds (Union[int, Tuple[int]]): length boundaries within which filtered sequences must be included
    - quality_threshold (int): reading quality value according to table phred+33, below which filtering will not be performed
    - output_filename (str): name of file where to save results without extension

    Return:

    - Dict[str, Tuple[str]]: a dictionary with filtered sequences

3. **`select_genes_from_gbk_to_fasta`**

Writes found CDS' in .fasta extension.

    Args:

    - input_gbk (str): input filename without .gbk extension
    - output_fasta (str): output filename without .fasta extension

    Return:

    - Dict[str, List[str]]: a dictionary with name of gene and its translation info

### Input of data

Functions `run_ultimate_protein_tools` and `run_dna_rna_tools` checks the input for mistakes and says if there are some of that, returning None for a specific sequence. 

Function `filter_fastq_seqs` and other functions, that work with files, takes as input the file name without extension!

During each run of the first 2 functions, the user is required to enter a **protein sequence / sequences** or **DNA / RNA sequence / sequences** that must be processed using the procedures listed above.

The program involves the analysis of protein sequences consisting of <u>**20 canonical amino acids**</u> and nucleic acids consisting of <u>**5 nucleotides (A, T, G, C, U)**</u>.

```python
ut.run_ultimate_protein_tools('get_protein_rnas_number', 'HJKASDKHSJAD')
'Some of your sequences have mistakes!'
None
```

### Examples

You can use the `EXAMPLE_FASTQ` constant for your tests!

```python
ut.filter_fastq_seqs('example_fastq')  # no parameteres = all sequences

ut.filter_fastq_seqs('example_fastq', gc_bounds=(50, 50), length_bounds=75, quality_threshold=32)  # if you want to indicate a specific sequence size or GC composition, then indicate it twice as an interval in parentheses

ut.run_dna_rna_tools('CGGGTCATTTTAAAGTGTGGCTAACTAGTAATATTCAAAGATGTCTTAC', 'reverse_complement')
'GTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG'

ut.run_dna_rna_tools('ATG', 'transcribe')
'AUG'

ut.run_dna_rna_tools('ATG', 'reverse')
'GTA'

ut.run_dna_rna_tools('AtG', 'complement')
'TaC'

ut.run_dna_rna_tools('ATG', 'aT', 'reverse')
['GTA', 'Ta']

ut.run_ultimate_protein_tools('get_fracture_of_aa', 'MMMfasdLA')
{'a': 0.1111, 'L': 0.1111, 'M': 0.3333, 'A': 0.1111, 's': 0.1111, 'd': 0.1111, 'f': 0.1111}

run_ultimate_protein_tools('count_aa', 'MMMfasdLA')
{'L': 1, 'M': 3, 'd': 1, 'f': 1, 'a': 1, 'A': 1, 's': 1}

ut.run_ultimate_protein_tools('is_protein_valid', 'MMMfasdLA', 'XXXDDDD')
'Some of your sequences have mistakes!'
[True, None]

bfp.convert_multiline_fasta_to_oneline('example_multiline_fasta')
"Output_fasta wasn't provided."
"<input_fasta>_oneline created!"

bfp.select_genes_from_gbk_to_fasta('example_bfp.bfp')
"Output_fasta wasn't provided."
"<input_gbk>.fasta created!"
```

### Contacts
Artem Vasilev (artem_vasilev_01@list.ru)
