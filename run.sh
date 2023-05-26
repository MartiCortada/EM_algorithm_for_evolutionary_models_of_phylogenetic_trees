#!/bin/bash

# Command 1
echo "running EM_a1_b1_500_random:"
python3 EM.py newick.txt random 500 100 EM_a1_b1_500_random

# Command 2
echo "running EM_a1_b1_1000_random:"
python3 EM.py newick.txt random 1000 100 EM_a1_b1_1000_random

# Command 3
echo "running EM_a1_b1_10000_random:"
python3 EM.py newick.txt random 10000 100 EM_a1_b1_10000_random

# Command 4
echo "running EM_a15_b15_500_random:"
python3 EM.py newick2.txt random 500 100 EM_a15_b15_500_random

# Command 5
echo "running EM_a15_b15_1000_random:"
python3 EM.py newick2.txt random 1000 100 EM_a15_b15_1000_random

# Command 6
echo "running EM_a15_b15_10000_random:"
python3 EM.py newick2.txt random 10000 100 EM_a15_b15_10000_random

# Command 7
echo "running EM_a30_b30_500_random:"
python3 EM.py newick3.txt random 500 100 EM_a30_b30_500_random

# Command 8
echo "running EM_a30_b30_1000_random:"
python3 EM.py newick3.txt random 1000 100 EM_a30_b30_1000_random

# Command 9
echo "running EM_a30_b30_10000_random:"
python3 EM.py newick3.txt random 10000 100 EM_a30_b30_10000_random

# Add more commands as needed
echo "All commands executed."