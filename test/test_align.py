# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    _, _, _, matrix = nw.align(seq1, seq2)

    expected_matrix = np.array([
        [0, -np.inf, -np.inf, -np.inf],
        [-np.inf, 5, -11, -13],
        [-np.inf, -12, 4, -8],
        [-np.inf, -12, -1, 5],
        [-np.inf, -14, -6, 4]
    ])
    
    np.testing.assert_array_equal(expected_matrix, matrix)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, aligned_seq1, aligned_seq2, _ = nw.align(seq3, seq4)

    assert score == 17, f"Expected score 17, got {score}"

    expected_seq1 = "MAVHQLIRRP"
    expected_seq2 = "M---QLIRHP"

    assert aligned_seq1 == expected_seq1, f"Expected {expected_seq1}, got {aligned_seq1}"
    assert aligned_seq2 == expected_seq2, f"Expected {expected_seq2}, got {aligned_seq2}"

if __name__ == "__main__":
    pytest.main()