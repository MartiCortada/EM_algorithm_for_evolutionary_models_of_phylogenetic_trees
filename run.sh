#!/bin/bash

# Command 1
echo "running EM_a75_b1_10000_random:"
python3 EM.py newick.txt random 10000 100 EM_a75_b1_10000_random

# Command 2
echo "running EM_a1_b75_500_random:"
python3 EM.py newick2.txt random 500 100 EM_a1_b75_500_random

# Command 3
echo "running EM_a1_b75_1000_random:"
python3 EM.py newick2.txt random 1000 100 EM_a1_b75_1000_random

# Command 4
echo "running EM_a1_b75_10000_random:"
python3 EM.py newick2.txt random 10000 100 EM_a1_b75_10000_random

# Add more commands as needed
echo "All commands executed."