"""
Algorithm maximum matching RNA secondary structures.
Problem: http://rosalind.info/problems/mmch/
Author: Sophie Hospel
Last update: 07-12-2020
"""


def get_file():
    """
    Ask user to enter filename which will be used for the calculation
    of the maximum matching.
    :return filename : Name of the Fasta file as String.
    """

    filename = input('Type filename: ')
    return filename


def read_file(filename):
    """
    Reads your Fasta file and splits the header from sequence. For each
    duo of a header and sequence the rest of the code will calculate and
    print the result. If there are newlines in the sequence, the sequence
    lines will be added to each other.
    :param filename: Name of the Fasta file as String.
    :return headers_and_sequences: Dict with the header as key and
    the sequence as value.
    """

    headers_and_sequences = {}

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line.split('\n')[0]
                sequence = ''
            elif line.startswith('A') or line.startswith('U') or line.startswith('G') or line.startswith('C'):
                sequence += line.split('\n')[0]
                headers_and_sequences[header] = sequence
            else:
                print("Can't find header in your file.")
                raise SystemExit

    return headers_and_sequences


def max_matching(sequence):
    """
    Calculates the maximum amount of matching possibilities. This is
    calculated by the count of appearance of each nucleotide. The nucleotides
    A and U are multiplied and G and C are multiplied. This output is multiplied
    with each other. This is the maximum amount of matching possibilities.
    :param sequence: Sequence found in Fasta file as String.
    :return result: Amount of possibilities A-U & G-C matches as Integer.
    """

    amount_a = 0
    amount_u = 0
    amount_g = 0
    amount_c = 0

    for nucleotide in sequence:
        if nucleotide.upper() == 'A':
            amount_a += 1
        elif nucleotide.upper() == 'U':
            amount_u += 1
        elif nucleotide.upper() == 'G':
            amount_g += 1
        elif nucleotide.upper() == 'C':
            amount_c += 1

    a_to_u = amount_a * amount_u
    g_to_c = amount_g * amount_c

    result = a_to_u * g_to_c

    return result


def print_result(header, sequence, result):
    """
    Prints the result that was calculated in the previous function.
    To add full information the header and sequence are added to the
    result String that is shown for the user in the console.
    :param header: Header found in Fasta file as String.
    :param sequence: Sequence found in Fasta file as String.
    :param result: Amount of possibilities A-U & G-C matches as Integer.
    """

    print('For header {0} with sequence {1} are {2} maximum matching results.'.format(header, sequence, result))


if __name__ == '__main__':
    filename = get_file()
    headers_and_sequences = read_file(filename)
    for header, sequence in headers_and_sequences.items():
        result = max_matching(sequence)
        print_result(header, sequence, result)