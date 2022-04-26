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

#define DEBUG false
#define debug_print(x) do {if (DEBUG)   std::cerr << x; } while (0)

#define NIL -1
#define ROOT_NODE 0
#define DEFAULT_TAG 0
#define ROOT_ID 1
#define BACKWARD_NODES 2
#define N_CHILDREN 2
#define EVEN_OFFSET_START 3
#define LAST_EDGE 4

unsigned log2base(unsigned val);

void mpi_error();

bool even(unsigned x);

template<typename T>
void print_single_arr(std::vector <T> arr, std::string name);

void prints_adjacent_representation(std::vector<unsigned> target_nodes,
                                    std::vector<int> next_edge, std::vector<unsigned> nodes_first_edge);

void load_adjacent_representation(int rank, unsigned *target_node, int *next_edge,
                                  unsigned tree_depth, unsigned n_edges);

int get_reversed_next_edge(unsigned rank, unsigned reversed_rank, int next_edge);

#endif //PREORDER_H
