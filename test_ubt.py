import os

import pytest

from bio_files_processor import convert_multiline_fasta_to_oneline, OpenFasta
from genscan_scripts import URL, ORGANISMS, EXON_CUTOFFS
from ultimate_tools import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    WrongAlphabetError,
)


@pytest.fixture
def tmp_file():
    file_path = "tmp.fasta"
    yield file_path
    if os.path.exists(file_path):
        os.remove(file_path)


class TestGenscan:
    def test_genscan_url(self):
        target = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"
        result = URL
        assert result == target

    def test_genscan_organisms(self):
        target = ["Vertebrate", "Arabidopsis", "Maize"]
        result = ORGANISMS
        assert result == target

    def test_genscan_exon_cutoffs(self):
        target = [1.00, 0.50, 0.25, 0.10, 0.05, 0.02, 0.01]
        result = EXON_CUTOFFS
        assert result == target


class TestSeqClasses:
    def test_dna_alphabet(self):
        seq = "ATGCU"
        with pytest.raises(WrongAlphabetError):
            DNASequence(seq).is_alphabet_correct()

    def test_rna_alphabet(self):
        seq = "ATGCU"
        with pytest.raises(WrongAlphabetError):
            RNASequence(seq).is_alphabet_correct()

    def test_protein_alphabet(self):
        seq = "FLSYZ"
        with pytest.raises(WrongAlphabetError):
            AminoAcidSequence(seq).is_alphabet_correct()


class TestFastaFiles:

    def test_oneline_fasta_exists(self, tmp_file):
        convert_multiline_fasta_to_oneline("data/example.fasta", tmp_file)
        assert os.path.exists(tmp_file)

    def test_oneline_fasta_content(self):
        with OpenFasta("data/example.fasta") as fasta:
            target_seqs = fasta.read_records()

        convert_multiline_fasta_to_oneline("data/example.fasta")
        with OpenFasta("data/example_oneline.fasta") as fasta:
            result_seqs = fasta.read_records()
        os.remove("data/example_oneline.fasta")
        assert target_seqs == result_seqs
