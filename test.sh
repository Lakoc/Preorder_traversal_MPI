#!/bin/bash

processes=12
input="ABCDEFG"

#preklad cpp zdrojaku
mpic++ -o preorder preorder.cpp

#spusteni
mpirun  -np ${processes} ./preorder ${input}

#uklid
rm -f preorder