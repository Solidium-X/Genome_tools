"""

File: Class - Sequence
Interpreter: Python 3.7
Author: Sebastian Pineda
E-mail: sebastian.pineda@hotmail.co.uk

Description:
    This class has various methods to help analyse DNA and RNA strands. It can be used to get information on separate
    genes of a genome. It can also be used to gain simple data from the strands such as base count and whether it is a
    valid DNA or RNA strand.

"""

from collections import Counter


class Sequence:

    # Private static class variables, do NOT alter!
    _is_valid_result = False
    _dna_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    _rna_complement = {'A': 'T', 'T': 'A', 'U': 'A', 'C': 'G', 'G': 'C'}
    _is_dna = False
    _is_rna = False
    _dna_viable_bases = ['A', 'T', 'C', 'G']
    _rna_viable_bases = ['A', 'T', 'C', 'G', 'U']

    def __init__(self, header, sequence, strand=None):
        """

        Args:
            header: Takes in a 'title header' for the DNA sequence.
            sequence: Takes in a DNA sequence as string, converts to a list.
                      Automatically converts all characters to capitals.
            strand: Can be set to 'DNA' or 'RNA', will affect how some class methods function.
                    Defaults to DNA if not specified.

        """
        self.header = header
        self.sequence = list(sequence.upper())
        self.strand = strand

        if strand in [None, 'dna', 'DNA']:
            self._is_dna = True
        elif strand in ['rna', 'RNA']:
            self._is_rna = True
        else:
            raise Exception("Unrecognised type: Specify either 'DNA' or 'RNA' strand.")

    def __str__(self):
        """

        this method is invoked to get a string representation of the object,
        e.g. if we want to print(..) it.

        """
        return '<' + self.header + '>\n' + self.sequence.__str__() + '\n'

    def __eq__(self, other):
        """

        Args:
            other: Takes in other.sequence.

        Returns:
            True: if other is of the same class AND the sequence is of equal value.
            False: otherwise.

        """
        if isinstance(other, self.__class__):
            return self.sequence == other.sequence
        else:
            return False

    def count_bases(self, style=None):
        """

        Method is altered if self._is_dna/rna is True. Returns error if neither are True.

        Args:
            self: Takes in self.sequence and looks for DNA bases.
            style: Can specify the style of the return as (i)ndividual, (t)otal or (b)oth. Defaults to total.

        Returns:
            Total and/or individual count of bases that are either 'A', 'T', 'C' or 'G' as an integer.

        """
        _count = Counter(self.sequence)

        if self._is_dna is True:
            if style in [None, 't', 'total']:
                return _count['A'] + _count['T'] + _count['C'] + _count['G']

            elif style in ['i', 'individual']:
                return _count

            elif style in ['b', 'both']:
                return _count, _count['A'] + _count['T'] + _count['C'] + _count['G']

        elif self._is_rna is True:
            if style in [None, 't', 'total']:
                return _count['A'] + _count['T'] + _count['C'] + _count['G'] + _count['U']

            elif style in ['i', 'individual']:
                return _count

            elif style in ['b', 'both']:
                return _count, _count['A'] + _count['T'] + _count['C'] + _count['G'] + _count['U']

        else:
            raise Exception('Error: class type (DNA or RNA) not specified.')

    def is_valid(self):
        """

        Method is altered if self._is_dna/rna is True. Returns error if neither are True.

        Args:
            self: Takes in self.sequence.

        Returns:
            Boolean value True if all elements in self.sequence contains only ['A', 'T', 'C', 'G' and/or 'U'] or
            boolean value False if otherwise.
            Saves results to static class variable (_is_valid_result).

        """
        if self._is_dna is True:
            if all(element in self._dna_viable_bases for element in self.sequence):
                self._is_valid_result = True
            return self._is_valid_result

        elif self._is_rna is True:
            if all(element in self._rna_viable_bases for element in self.sequence):
                self._is_valid_result = True
            return self._is_valid_result

    def dna_complement(self):
        """

        Method is altered if self._is_dna/rna is True. Returns error if neither are True.

        Args:
            self: Takes in self.sequence DNA strand.

        Returns:
            Complement strand strand as list (3' - 5').

        """
        if self._is_dna is True:
            return list("".join([self._dna_complement[base] for base in self.sequence]))

        elif self._is_rna is True:
            return list("".join([self._rna_complement[base] for base in self.sequence]))

        else:
            raise Exception('Error: class type (DNA or RNA) not specified.')

    def dna_reverse_complement(self):
        """

        Method is altered if self._is_dna/rna is True. Returns error if neither are True.

        Args:
            self: takes in self.sequence DNA strand.

        Returns:
            Reverse complement strand as list (5' - 3').

        """
        if self._is_dna is True:
            return list("".join([self._dna_complement[base] for base in reversed(self.sequence)]))

        elif self._is_rna is True:
            return list("".join([self._rna_complement[base] for base in reversed(self.sequence)]))

        else:
            raise Exception('Error: class type (DNA or RNA) not specified.')

    def find_mutation(self, other):
        """

        Args:
            self: takes in self.sequence DNA strand.
            other: Takes in other.sequence and compares to self.sequence

        Returns:
            length exception: if both sequences are not of the same length, raises an exception error.
            index: if both sequence lengths match and if there is a base mutation, returns the first occurring
            mutation index as integer.
            -1: if both sequence lengths AND bases match completely, returns integer.
            error exception: this exception should not raise if the code is running correctly. Placed for debugging.

        """
        if len(self.sequence) != len(other.sequence):
            raise Exception('Cannot compare sequences of different lengths.')

        elif len(self.sequence) == len(other.sequence):
            print('sequence lengths match')

            for index, (i, j) in enumerate(zip(self.sequence, other.sequence)):
                if i == 'A' and j != 'A' \
                        or i == 'T' and j != 'T' \
                        or i == 'C' and j != 'C' \
                        or i == 'G' and j != 'G':
                    return index

            else:
                return -1

        else:
            raise Exception("Error: oops that shouldn't have happened!")

    def split_genes(self):
        """

        Args:
            Takes in self.sequence.

        Returns:
            self: if no full-stop signal is detected in the sequence, it will return self.
            gene_list: if full-stop signal is detected, it will return list of class instances for per gene.

        """
        full_stop = 'AAAAAAAAAATTTTTTTTTT'
        seq = ''.join(self.sequence)

        if full_stop not in seq:
            print('This genome does not contain full-stop signals')
            return self

        elif full_stop in seq:
            gene_list = []
            gene_count = 1
            genes = seq.split(full_stop)
            for gene in genes:
                cls_instance = Sequence(f'Gene {gene_count}  from {self.header}', gene)
                gene_list.append(cls_instance)
                gene_count += 1
            # print(f'Total number of genes for this genome = {gene_count-1}')
            return gene_list

    def count_swaps(self, other):
        """

        Args:
            self: Takes in self.sequence.
            other: Takes in other.sequence for comparison to self.sequence.

        Returns:
            length exception: if length of self.sequence and other.sequence are different
            count: if self and other length matches, returns number of swap mutations as integer.
            error exception: this exception should not raise if the code is running correctly. Placed for debugging.

        """
        if len(self.sequence) != len(other.sequence):
            raise Exception('Cannot compare sequences of different lengths')

        elif len(self.sequence) == len(other.sequence):
            count = 0

            for index, (i, j) in enumerate(zip(self.sequence, other.sequence)):
                if i == 'A' and j != 'A' \
                        or i == 'T' and j != 'T' \
                        or i == 'C' and j != 'C' \
                        or i == 'G' and j != 'G':
                    count += 1
            return count

        else:
            raise Exception("Error: oops that shouldn't have happened!")

    def base_count_ingenes(self):
        """

        Args:
            self: Takes in self.sequence (typically a genome). Subjects self.sequence to .split_gene() class method.

        Returns:
            List containing base counts as integers for each gene of self.sequence.

        """
        gene_list = self.split_genes()
        iteration = 0
        gene_length_list = []
        for _ in gene_list:
            gene_length = gene_list[iteration].count_bases()
            gene_length_list.append(gene_length)
            iteration += 1
        return gene_length_list

    def compare_swaps_ingenes(self, other):
        """

        Args:
            self: Takes self.sequence and applies .split_genes() class method to it.
            other: Takes other.sequence and applies .split_genes() class method to it. Compared to self.

        Returns:
            mutation_list: A list of integers containing swap mutations for per gene, when self compared to other.

        """
        self_gene_list = self.split_genes()
        other_gene_list = other.split_genes()
        mutation_list = []
        index = 0
        for _ in self.split_genes():
            mutations = self_gene_list[index].count_swaps(other_gene_list[index])
            mutation_list.append(mutations)
            index += 1

        return mutation_list


def seq_read_file(filename):
    """

    Args:
        filename: the relative path of the file that is to be opened and read.

    Returns:
        An instance of the class Sequence using the data contained within filename.

    """
    with open(filename, 'r', encoding='ASCII') as dnafile:
        file_contents = dnafile.read()
        content_split = file_contents.split('\n')
        header = content_split[0]
        genome = content_split[1]

        return Sequence(header, genome)
