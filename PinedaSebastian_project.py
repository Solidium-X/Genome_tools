"""

File: Test script for Class - Sequence
Interpreter: Python 3.7
Author: Sebastian Pineda
E-mail: sebastian.pineda@hotmail.co.uk

Description:
    This is a test script to demonstrate the various class methods and functions of the class - Sequence. Sections
    within the terminal are separated by headers to make it easier to read.

"""

from Final_project.PinedaSebastian_project_sequence import Sequence, seq_read_file
import matplotlib.pyplot as plt


# Task 1, 2 & 11 set up (NOTE: this also takes care of task 11)
print('#'*30, 'TASK 1, 2 & 11', '#'*30)
test_lower_sequence = 'gtcaTGca'   # lower case letters to see if it converts to upper case

# Task 1: create instance of class.
test_lower = Sequence(header='Lower header test', sequence=test_lower_sequence)

# Task 11: test to see if lower case sequence is converted to capitals when called through class.
print(test_lower.header)    # see if it prints the header
print(test_lower.sequence)  # see if it prints the sequence as list
print(test_lower)           # print the whole of the sequence with header on top as specified in __str__.

# Task 2 test using default settings (same as total style)
print('default: ', test_lower.count_bases())
assert test_lower.count_bases() == 8

# Task 2 test using individual style.
print('individual test: ', test_lower.count_bases(style='i'))
print('individual test: ', test_lower.count_bases(style='individual'))

# Task 2 test using total style.
print('total test: ', test_lower.count_bases(style='t'))
print('total test: ', test_lower.count_bases(style='total'))

# Task 2 test using both styles.
print('both test: ', test_lower.count_bases(style='b'))
print('both test: ', test_lower.count_bases(style='both'))


# Task 3 set up.
print('#'*30, 'TASK 3', '#'*30)
ts3_incorrect_sequence = 'gtcaKPca'     # contains 'K' and 'P'
ts3_test_incorrect = Sequence(header='Incorrect test', sequence=ts3_incorrect_sequence)

# Task 3 test using correct sequence from task 2. Should be True.
assert test_lower.is_valid() is True
# Task 3 test using incorrect sequence. Should be False.
assert ts3_test_incorrect.is_valid() is False
# Task 3 test to see if result of method is saved and can be recalled. Should be True.
assert test_lower._is_valid_result is True
print('See asserts within script')

# Task 4 set up. Both variables have same value
print('#'*30, 'TASK 4', '#'*30)
ts4_seq_1 = 'ATCGATCG'
ts4_seq_2 = 'ATCGATCG'

# Task 4 set up. Both variables used to create different instances of class.
ts4_seq_x = Sequence('task 4 sequence 1', ts4_seq_1)
ts4_seq_y = Sequence('task 4 sequence 2', ts4_seq_2)

# Task 4 test
assert ts4_seq_y == ts4_seq_x               # Same class, different variables but same value.
assert ts4_seq_x == ts4_seq_x               # Comparing to self.
assert ts4_seq_1 != ts4_seq_x.sequence      # 1 class and 1 variable with same value.
assert ts4_seq_1 != ts4_seq_y.sequence      # Same as above.
print('See asserts within script')

# Task 5 set up.
print('#'*30, 'TASK 5', '#'*30)

ts5_seq = 'GGATCAAAGCACTTCGACTC'
ts5_seq_complement = 'CCTAGTTTCGTGAAGCTGAG'
ts5_seq_reverse = 'GAGTCGAAGTGCTTTGATCC'

ts5_cls = Sequence('task 5 sequence', ts5_seq)
ts5_cls_complement = Sequence('task 5 complement', ts5_seq_complement)
ts5_cls_reverse = Sequence('task 5 reverse', ts5_seq_reverse)

# Task 5 test
print(ts5_seq)                                              # print to see sequence is formed correctly
print(ts5_cls)                                              # print whole class.
print('Complement: ', ts5_cls.dna_complement())             # see if complement is correct.
print('Reverse: ', ts5_cls.dna_reverse_complement())        # see if complement is correct when reversed.

# Task 5 comparison tests
print('Compare complement: ', ts5_cls.dna_complement())
print('Compare complement: ', ts5_cls_complement.sequence)
assert ts5_cls_complement.sequence == ts5_cls.dna_complement()      # Check to see if calculated complement is correct.
assert ts5_cls_complement.sequence != ts5_cls_reverse.sequence      # Same as above but complement reversed.


# Task 6 set up.
print('#'*30, 'TASK 6', '#'*30)
ts6_seq1_match = ('A' * 10000)
ts6_seq2_match = ('A' * 10000)
ts6_seq1_mutation = ('A' * 10000) + 'A'
ts6_seq2_mutation = ('A' * 10000) + 'C'

# Task 6 classing of variables.
ts6_cls1 = Sequence('task6 sequence 1', ts6_seq1_match)
ts6_cls2 = Sequence('task6 sequence 2', ts6_seq2_match)
ts6_cls1_mut = Sequence('task6 sequence 1 mutation', ts6_seq1_mutation)
ts6_cls2_mut = Sequence('task6 sequence 2 mutation', ts6_seq2_mutation)

# Task 6 match testing.
compare = ts6_cls1.find_mutation(ts6_cls2)
# Should return -1
print(compare)
assert compare == -1

# Task 6 mutation testing.
compare_mis = ts6_cls1_mut.find_mutation(ts6_cls2_mut)
# Should return 0-index of first mismatching base-pair (in this case 10000).
print(compare_mis)
assert compare_mis == 10000


# Task 7 set up
print('#'*30, 'TASK 7', '#'*30)

# read file and initiate new instance of Sequence using file.
ts7_test = seq_read_file('genome_01.dat')

# print total bases of sequence using class function.
print('Total number of bases for genome_01.dat =', ts7_test.count_bases())


# Task 8 set up
print('#'*30, 'TASK 8', '#'*30)
non_stop_sequence = 'CCGATCGAAA'
example_sequence = 'CCGATCGAAAAAAAAAAATTTTTTTTTT'

ts8_cls_1 = Sequence('task 8 no stops', non_stop_sequence)
ts8_cls_2 = Sequence('task 8 example', example_sequence)

# Test if method handles genome with no full-stop signals. Should return the class instance.
ts8_list_1 = ts8_cls_1.split_genes()
print(ts8_list_1)

# Test if method can handle the example given in instructions and return the first gene.
ts8_list_2 = ts8_cls_2.split_genes()
print(ts8_list_2[0])
# Test if class instance within list can be called using class method. Also returns length of first gene.
print('Length of gene 1 = ', ts8_list_2[0].count_bases(), '\n')
assert ts8_list_2[0].count_bases() == 8

# Test to see if method can handle file read, separating genes and assigning the class instances to a list.
ts8_list_3 = ts7_test.split_genes()
print('File open test here: \n', ts8_list_3[0])
# Test to see if class instances within list can be called using class methods. Also returns length of first gene.
print('Length of gene 1 = ', ts8_list_3[0].count_bases())


# Task 9 set up
print('#'*30, 'TASK 9', '#'*30)
ts9_test = seq_read_file('genome_01.dat')

# This class method auto splits the class instance sequence into genes and returns the genes into a list.
ts9_genelength_list = ts9_test.base_count_ingenes()

# Check to see if gene lengths were extracted correctly
print(ts9_genelength_list)

# Plotting of gene length data to histograms
x = ts9_genelength_list[:-1]      # excluding the last 'gene' produces a cleaner histogram.
plt.hist(x, bins=25)
plt.ylabel('Occurrence')
plt.xlabel('Gene length (base pairs)')
plt.title('Histogram for occurrence of gene lengths of genes in genome_02.dat')
plt.savefig('Task 9 histogram.png')
plt.show()


# Task 10 setup
print('#'*30, 'TASK 10', '#'*30)
ts10_cls_1 = seq_read_file('genome_01.dat')
ts10_cls_2 = seq_read_file('genome_02.dat')

# Test to see if it can count mutation swaps between the two genomes (takes a while to calculate)
print('Mutation swaps between genome_01 and genome_02: ', ts10_cls_1.count_swaps(ts10_cls_2))

# Get base counts per genes.
ts10_list1_basecounts = ts10_cls_1.base_count_ingenes()
ts10_list2_basecounts = ts10_cls_2.base_count_ingenes()
# Check they give the same value, meaning each gene is the same length.
assert ts10_list1_basecounts == ts10_list2_basecounts

# Instances of Sequence automatically have .split_genes applied to self and then have the gene lists compared to
# argument (other.split_genes) to provide the number of mutations per gene.
mutation_list = ts10_cls_1.compare_swaps_ingenes(ts10_cls_2)

# Check code is working correctly
print(mutation_list)
assert sum(mutation_list) == 272

# Plotting of scatter-plot
x = ts10_list2_basecounts[:-1]      # exclude last 'gene' to make graph cleaner
y = mutation_list[:-1]              # exclude last 'gene' to make graph cleaner

plt.scatter(x, y, alpha='0.5')
plt.xlabel('Gene length (base-pairs)')
plt.ylabel('Number of swap mutations')
plt.title('Scatter plot of number of mutations against gene length')
plt.savefig('Task 10 scatterplot.png')
plt.show()


# Task 13 setup
print('#'*30, 'TASK 13', '#'*30)
ts13_rna = 'ATCGATCGUATCGU'                 # Should return 14
ts13_rna_incorrect = 'ATCGATCGUATCGUPP'     # Contains 'P' should not return True when is_valid() is called.
ts13_dna = 'ATCGATCGATCG'                   # Should return 12
ts13_dna_incorrect = 'ATCGATCGATCGPP'       # Contains 'P' should not return True when is_valid() is called.

ts13_rna_cls = Sequence('Task 13 RNA', ts13_rna, strand='RNA')
ts13_rna_cls_incorrect = Sequence('Task 13 RNA', ts13_rna_incorrect, strand='RNA')
ts13_dna_cls = Sequence('Task 13 DNA', ts13_dna, strand='DNA')
ts13_dna_cls_incorrect = Sequence('Task 13 DNA', ts13_dna_incorrect, strand='DNA')

# Should pass without error to validate that .isvalid() is working correctly.
assert ts13_rna_cls.is_valid() is True
assert ts13_dna_cls.is_valid() is True
assert ts13_rna_cls_incorrect.is_valid() is False
assert ts13_dna_cls_incorrect.is_valid() is False
print('See asserts within script')
