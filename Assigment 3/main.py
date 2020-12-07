"""
Implementing Needleman-Wunsch global sequence alignment algorithm.
Author: Sophie Hospel
Last update: 07-12-2020
"""

import nwalign3 as nw


def get_alignment_strings():
    """
    Ask user to enter alignment Strings which will be used for the
    global alignment.
    :return alignment1: First String for alignment.
    :return aligment2: Second String for alignment.
    """

    alignment1 = input('Enter first word you want to align: ')
    alignment2 = input('Enter first word you want to align: ')

    return alignment1, alignment2


def align_and_print(alignment1, alignment2):
    """
    Aligns the two given words by the user. Matrix is set to PAM250.
    Prints the result of the alignment.
    :param alignment1: First String for alignment.
    :param alignment2: Second String for alignment.
    :return:
    """

    global_align = nw.global_align(alignment1, alignment2, matrix='PAM250')
    print('Global aligment', global_align)


if __name__ == '__main__':
    alignment1, alignment2 = get_alignment_strings()
    align_and_print(alignment1, alignment2)
