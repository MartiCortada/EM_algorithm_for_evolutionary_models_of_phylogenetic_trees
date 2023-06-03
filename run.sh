#!/bin/bash

# Command 1
echo "running EM_a10_b10_500_random:"
python3 EM_com.py newick.txt 500 100

# Command 2
echo "running EM_a10_b10_1000_random:"
python3 EM_com.py newick.txt 1000 100

# Command 3
echo "running EM_a10_b10_10000_random:"
python3 EM_com.py newick.txt 10000 100

# Add more commands as needed
echo "All commands executed."