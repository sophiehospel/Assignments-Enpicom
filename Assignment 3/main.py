"""
Implementing the Needleman-Wunsch algorithm.
Author: Sophie Hospel
Last update: 08-12-2020
"""

gap_penalty = -1
match_award = 1
mismatch_penalty = -1


def get_sequences():
    """
    Ask user to enter sequences which will be used for the alignment.
    :return seq1: Sequence one as String.
    :return seq2: Sequence two as String.
    """

    seq1 = input('Type first sequence: ')
    seq2 = input('Type second sequence: ')

    return seq1, seq2


def create_zero_matrix(length_seq2, length_seq1):
    """
    Creates a matrix for the frame, filled with zeros. Add sublists
    in the amount of the length sequence1 + 1. Add zeros
    to the sublists in the amount of the length sequence2 + 1.
    :param length_seq2: Length of sequence 2 + 1.
    :param length_seq1: Length of sequence 1 + 1.
    :return zero_matrix: List with sublists in the amount of sequence 2,
    in the sublists zeros in the amount of sequence 1.
    """

    zero_matrix = []

    for x in range(length_seq2):
        zero_matrix.append([])
        for y in range(length_seq1):
            zero_matrix[-1].append(0)

    return zero_matrix


def create_score_matrix(length_seq1, length_seq2, zeros_matrix):
    """
    Calculates the scores for in the matrix. First fill out the horizontal
    column from 0 to -(length sequence 2). Then fill out the vertical column
    from 0 to -(length sequence 1). After those steps, fill out the other
    values in the zeros_matrix matrix. Calculate the zeros_matrix by checking
    the top, left, and diagonal cells. Record the maximum zeros_matrix from
    the three possible scores calculated (match, delete and insert).
    :param length_seq1: Length of sequence 1.
    :param length_seq2: Length of sequence 2.
    :param zeros_matrix: List with sublists in the amount of sequence 2,
    in the sublists zeros in the amount of sequence 1.
    :return zeros_matrix: List with sublists in the amount of sequence 2,
    in the sublists the best possible alignment score of the sequences in the
    amount of sequence 1.
    """

    for i in range(0, length_seq2 + 1):
        zeros_matrix[i][0] = gap_penalty * i

    for j in range(0, length_seq1 + 1):
        zeros_matrix[0][j] = gap_penalty * j

    for i in range(1, length_seq2 + 1):
        for j in range(1, length_seq1 + 1):
            match = zeros_matrix[i - 1][j - 1] + match_score(seq1[j - 1], seq2[i - 1])
            delete = zeros_matrix[i - 1][j] + gap_penalty
            insert = zeros_matrix[i][j - 1] + gap_penalty

            zeros_matrix[i][j] = max(match, delete, insert)

    return zeros_matrix


def match_score(nucleotide1, nucleotide2):
    """
    A function for determining the score between any two bases in alignment.
    :param nucleotide1: Nucleotide of sequence 1.
    :param nucleotide2: Nucleotide of sequence 2.
    :return gap_penalty: Global variable set to a gap penalty score.
    :return mismatch_penalty: Global variable set to a mismatch penalty score.
    """

    if nucleotide1 == nucleotide2:
        return match_award
    elif nucleotide1 == '-' or nucleotide2 == '-':
        return gap_penalty
    else:
        return mismatch_penalty


def needleman_wunsch(seq1, seq2):
    """
    Implementing the Needleman-Wunsch algorithm. Gets the frame for the matrix with
    function create_zero_matrix(). Fills the matrix with scores with function
    create_score_matrix(). Calculates the best possible alignment by the scores.
    The lines are in the calculation reverse for the characters in each sequence.
    Set forward to show as result and save.
    :param seq1: Sequence one as String.
    :param seq2: Sequence two as String.
    :return alignment1: first best possible alignment of sequence 1.
    :return alignment2: first best possible alignment of sequence 2.
    """

    length_seq1 = len(seq1)
    length_seq2 = len(seq2)

    zeros_matrix = create_zero_matrix(length_seq2 + 1, length_seq1 + 1)
    score_matrix = create_score_matrix(length_seq1, length_seq2, zeros_matrix)

    alignment1 = ''
    alignment2 = ''

    i = length_seq2
    j = length_seq1

    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i][j - 1]
        score_left = score_matrix[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[j - 1], seq2[i - 1]):
            alignment1 += seq1[j - 1]
            alignment2 += seq2[i - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            alignment1 += seq1[j - 1]
            alignment2 += '-'
            j -= 1
        elif score_current == score_left + gap_penalty:
            alignment1 += '-'
            alignment2 += seq2[i - 1]
            i -= 1

    while j > 0:
        alignment1 += seq1[j - 1]
        alignment2 += '-'
        j -= 1
    while i > 0:
        alignment1 += '-'
        alignment2 += seq2[i - 1]
        i -= 1

    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]

    return alignment1, alignment2


def print_result(seq1, seq2, alignment1, alignment2):
    """
    Formats the results from the alignment. To add full information
    the sequences with the algorithm name are printed. Then printed
    alignment1 and alignment2.
    :param seq1: Sequence one as String.
    :param seq2: Sequence two as String.
    :param alignment1: first best possible alignment of sequence 1.
    :param alignment2: first best possible alignment of sequence 2.
    """

    print()
    print('Your Needleman-Wunsch alignment of sequences {0} and {1}.'.format(seq1, seq2))
    print(alignment1)
    print(alignment2)


if __name__ == '__main__':
    seq1, seq2 = get_sequences()
    alignment1, alignment2 = needleman_wunsch(seq1, seq2)
    print_result(seq1, seq2, alignment1, alignment2)
