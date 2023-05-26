#!/bin/bash

# Command 1
echo "running EM_a30_b30_1000_random:"
python3 EM.py newick3.txt random 1000 100 EM_a30_b30_1000_random

# Command 2
echo "running EM_a30_b30_10000_random:"
python3 EM.py newick3.txt random 10000 100 EM_a30_b30_10000_random

# Command 1
echo "running EM_a45_b45_500_random:"
python3 EM.py newick.txt random 500 100 EM_a45_b45_500_random

# Command 1
echo "running EM_a45_b45_1000_random:"
python3 EM.py newick.txt random 1000 100 EM_a45_b45_1000_random

# Command 1
echo "running EM_a45_b45_10000_random:"
python3 EM.py newick.txt random 10000 100 EM_a45_b45_10000_random

# Add more commands as needed
echo "All commands executed."