#!/bin/bash

# Command 1
echo "running EM_a75_b75_1000_random:"
python3 EM.py newick2.txt random 1000 100 EM_a75_b75_1000_random

# Command 1
echo "running EM_a60_b60_1000_random:"
python3 EM.py newick.txt random 1000 100 EM_a60_b60_1000_random

# Add more commands as needed
echo "All commands executed."