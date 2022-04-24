#!/bin/bash

# get tree
input=$1

# count input len
input_len=${#input}

# get number of processes
processes=$(($input_len * 2 - 2))

# compile source files
mpic++ -o preorder preorder.cpp

# run program
mpirun -oversubscribe -np ${processes} ./preorder ${input}

# clean up
rm -f preorder