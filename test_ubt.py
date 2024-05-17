import pytest

from genscan_scripts import URL, ORGANISMS, EXON_CUTOFFS
from ultimate_tools import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    WrongAlphabetError,
)


def test_genscan_url():
    target = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"
    result = URL
    assert result == target


def test_genscan_organisms():
    target = ["Vertebrate", "Arabidopsis", "Maize"]
    result = ORGANISMS
    assert result == target


def test_genscan_exon_cutoffs():
    target = [1.00, 0.50, 0.25, 0.10, 0.05, 0.02, 0.01]
    result = EXON_CUTOFFS
    assert result == target


def test_dna_alphabet():
    seq = "ATGCU"
    with pytest.raises(WrongAlphabetError):
        DNASequence(seq).is_alphabet_correct()


def test_rna_alphabet():
    seq = "ATGCU"
    with pytest.raises(WrongAlphabetError):
        RNASequence(seq).is_alphabet_correct()


def test_protein_alphabet():
    seq = "FLSYZ"
    with pytest.raises(WrongAlphabetError):
        AminoAcidSequence(seq).is_alphabet_correct()
