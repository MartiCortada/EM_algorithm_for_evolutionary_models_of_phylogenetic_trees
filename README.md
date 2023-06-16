# Expectation-Maximization for evolutionary models of phylogenetic trees

In order to test the **EM.py** implementation of Expectation-Maximization, execute the following command-line:

```python3 EM.py <newick_file> <init_method> <length_alignment> <repetitions> <name_test>```

where ```init_mode``` can be either **random** or **spectral** (which follows Chang's spectral theorem).

For example, 

*python3 EM.py tree.txt random 1000 100 test_tree*

runs the EM over the Newick tree in text format of tree.txt, with alignments of length 1000, doing 100 repetitions and storing the results in a directory named test_tree.
