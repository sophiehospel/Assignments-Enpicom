"""
Algorithm to compare DNA sequences with point mutations.
Problem: http://rosalind.info/problems/corr/
Author: Sophie Hospel
Last update: 07-12-2020
"""


def get_file():
    """
    Ask user to enter filename which will be used for comparing
    the sequences.
    :return filename : Name of the Fasta file as String.
    """

    filename = input('Type filename: ')
    return filename


def read_file(filename):
    """
    Reads your Fasta file and splits the header from sequence.
    If there are newlines in the sequence, the sequence lines
    will be added to each other.
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
            elif line.startswith('A') or line.startswith('T') or line.startswith('G') or line.startswith('C'):
                sequence += line.split('\n')[0]
                headers_and_sequences[header] = sequence

    return headers_and_sequences


def compare_sequences(headers_and_sequences):
    """
    Compares sequences from your Fasta file. The same sequences will be saved.
    Even as the sequences with maximal 1 difference in the sequence.
    :param headers_and_sequences: Dict with the header as key and
    the sequence as value.
    :return result: List with headers and sequences which appear twice in
    the Fasta file, or appear with one point mutation.
    """

    compare_set = headers_and_sequences
    result = []

    for header1, seq1 in headers_and_sequences.items():
        for header2, seq2 in compare_set.items():
            if header1 != header2:
                if seq1 != seq2:
                    differences = 0
                    for i in range(0, len(seq1)):
                        if seq1[i] != seq2[i]:
                            differences += 1
                    if differences == 1:
                        if any(header1 in s for s in result) and any(header2 in s for s in result):
                            continue
                        else:
                            result.append(' '.join(['(', header1, header2, ')', seq1, '->', seq2]))
                else:
                    if any(header1 in s for s in result) and any(header2 in s for s in result):
                        continue
                    else:
                        result.append(' '.join(['(', header1, header2, ')', seq1, '->', seq2]))
            else:
                continue

    return result


def get_reverse_complement(headers_and_sequences):
    """
    Get reverse complements from the sequences in your Fasta
    file for comparing in further code.
    :param headers_and_sequences: Dict with the header as key and
    the sequence as value.
    :return rc_seqs: Dict with the header as key and
    the reverse complement sequence as value.
    """

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rc_seqs = {}

    for header, seq in headers_and_sequences.items():
        bases = list(seq)
        bases = reversed([complement.get(base, base) for base in bases])
        bases = ''.join(bases)
        rc_seqs[header] = bases

    return rc_seqs


def check_reverse_complement(headers_and_sequences, rc_seqs, result):
    """
    Compares sequences from your Fasta file with the reverse complement
    sequence from the Fasta file. The same sequences will be saved.
    Even as the sequences with maximal 1 difference in the sequence.
    :param headers_and_sequences: Dict with the header as key and
    the sequence as value.
    :param rc_seqs: Dict with the header as key and
    the reverse complement sequence as value.
    :return result: List with headers and sequences which appear
    twice in the Fasta file, or appear with one point mutation.
    """

    for header1, seq1 in headers_and_sequences.items():
        for header2, seq2 in rc_seqs.items():
            if header1 != header2:
                if seq1 != seq2:
                    differences = 0
                    for i in range(0, len(seq1)):
                        if seq1[i] != seq2[i]:
                            differences += 1
                    if differences == 1:
                        if any(header1 in s for s in result) and any(header2 in s for s in result):
                            continue
                        else:
                            result.append(' '.join(['(', header1, header2, ')', seq1, '->', seq2]))
            else:
                continue

    return result


def print_result(result):
    """
    Formats the results as asked. To add full information the headers are added to the
    result String that is shown for the user in the console.
    :param result: List with headers and sequences which appear
    twice (with reverse checked) in the Fasta file, or appear with
    maximal one point mutation.
    """

    for res in result:
        print(res)


if __name__ == '__main__':
    filename = get_file()
    headers_and_sequences = read_file(filename)
    result = compare_sequences(headers_and_sequences)
    rc_seqs = get_reverse_complement(headers_and_sequences)
    result = check_reverse_complement(headers_and_sequences, rc_seqs, result)
    print_result(result)
