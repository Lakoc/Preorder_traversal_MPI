#ifndef PREORDER_H
#define PREORDER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <string>
#include <vector>

#define DEBUG true
//#define debug_print(fmt, ...) do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
//#define debug_print_(fmt) do { if (DEBUG) fprintf(stderr, fmt); } while (0)
#define debug_print(x) do {if (DEBUG)   std::cerr << x; } while (0)

#define ROOT_NODE 1
#define BACKWARD_NODES 2
#define EVEN_OFFSET_START 3
#define N_CHILDREN 2
#define NILL -1

unsigned log2base(unsigned val);
void mpi_error();
bool even(unsigned x);

#endif //PREORDER_H
