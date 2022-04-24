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
#define debug_print(x) do {if (DEBUG)   std::cerr << x; } while (0)

#define ROOT_NODE 0
#define ROOT_ID 1
#define LAST_EDGE 4
#define BACKWARD_NODES 2
#define EVEN_OFFSET_START 3
#define N_CHILDREN 2
#define NIL -1
#define DEFAULT_TAG 0
#define REQUEST_TAG 1
#define NO_REQUEST -1

unsigned log2base(unsigned val);

void mpi_error();

bool even(unsigned x);

template<typename T>
void print_single_arr(std::vector<T> arr, std::string name);

void prints_adjacent_representation(std::vector<unsigned> target_nodes,
                                    std::vector<int> next_edge, std::vector<unsigned> nodes_first_edge);

void load_adjacent_representation(std::vector<unsigned> &target_nodes,
                                  std::vector<int> &next_edge, std::vector<unsigned> &nodes_first_edge,
                                  unsigned tree_depth, unsigned n_nodes, unsigned n_edges);

#endif //PREORDER_H
