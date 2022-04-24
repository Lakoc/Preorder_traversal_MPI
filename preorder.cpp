#include "preorder.h"

/*
 * Function:  mpi_error
 * --------------------
 * Exit program with failure due to MPI runtime err.
 *
 */
void mpi_error() {
    printf("MPI library operation failed.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

/*
 * Function:  log2base
 * --------------------
 * Calculate logarithm of base 2 of unsigned number
 *
 *  val: value
 *
 *  returns: unsigned number equal log2 val
 */
unsigned log2base(unsigned val) {
    unsigned result = 0;
    while (val >>= 1) ++result;
    return result;
}


/*
 * Function:  even
 * --------------------
 * Check whatever number is even
 *
 *  val: value
 *
 *  returns: bool
 */
bool even(unsigned x) {
    if (x % 2 == 0) {
        return true;
    } else {
        return false;
    }
}

/*
 * Function:  prints_adjacent_representation
 * --------------------
 * Print representation of adjacent matrix
 *
 *  target_nodes: array of edge target nodes
 *
 *  next_edge: array of `pointers` to the next edge, contains NIL if next edge does not exist
 *
 *  nodes_first_edge: array of first edges reachable from node
 */
void prints_adjacent_representation(std::vector<unsigned> target_nodes,
                                    std::vector<int> next_edge, std::vector<unsigned> nodes_first_edge) {
    print_single_arr(target_nodes, "target_nodes");
    print_single_arr(next_edge, "next_edge");
    print_single_arr(nodes_first_edge, "nodes_first_edge");
}

/*
 * Function:  print_single_arr
 * --------------------
 * Print single arr
 *
 *  arr: arr to be printed
 *
 *  name: name of array
 *
 */
template<typename T>
void print_single_arr(std::vector <T> arr, std::string name) {
    debug_print(name << ": ");
    for (T item: arr) {
        debug_print(item << " ");
    }
    debug_print(std::endl);
}

/*
 * Function:  load_adjacent_representation
 * --------------------
 * Create adjacent matrix in form of multiple arrays, fill up following arrays.
 *
 *  target_nodes: array of edge target nodes
 *
 *  next_edge: array of `pointers` to the next edge, contains NIL if next edge does not exist
 *
 *  nodes_first_edge: array of first edges reachable from node
 *
 *  tree_depth: depth of tree - 1
 *
 *  n_nodes: number of nodes in tree
 *
 *  n_edges: number of edges in tree
 */
void load_adjacent_representation(std::vector<unsigned> &target_nodes,
                                  std::vector<int> &next_edge, std::vector<unsigned> &nodes_first_edge,
                                  unsigned tree_depth, unsigned n_nodes, unsigned n_edges) {
    unsigned first_forward = ROOT_ID + 1;
    unsigned backward_start = ROOT_ID * N_CHILDREN;
    unsigned even_offset = EVEN_OFFSET_START;
    for (unsigned i = 0; i < n_edges; ++i) {
        if (even(i)) {
            target_nodes.push_back(first_forward);
            first_forward += 1;
            if (even(i / 2) && i + N_CHILDREN + 1 <= n_edges) {
                next_edge.push_back(i + N_CHILDREN + 1);
            } else {
                next_edge.push_back(NIL);
            }
        } else {
            target_nodes.push_back(backward_start / BACKWARD_NODES);
            backward_start += 1;
            if (i + even_offset + 1 <= n_edges) {
                if (tree_depth <= 1) {
                    next_edge.push_back(NIL);
                } else {
                    next_edge.push_back(i + even_offset + 1);
                }
                even_offset += N_CHILDREN;
            } else {
                next_edge.push_back(NIL);
            }
        }
    }

    nodes_first_edge.push_back(ROOT_ID);
    for (unsigned i = 1; i < n_nodes; ++i) {
        nodes_first_edge.push_back(i * N_CHILDREN);
    }
}

int main(int argc, char **argv) {
    if (MPI_Init(&argc, &argv)) { mpi_error(); }
    // Load input tree
    std::string input_tree = argv[1]; // Name of the current exec program

    // Calculate number of edges, nodes and tree depth
    unsigned n_nodes = input_tree.length();
    unsigned n_edges = (2 * n_nodes) - 2;
    unsigned tree_depth = ceil(log2(n_nodes));

    // Auxiliary arrays to keep adjacency matrix in linear complexity
    std::vector<unsigned> target_nodes;
    std::vector<int> next_edge;
    std::vector<unsigned> nodes_first_edge;

    // Auxiliary variables to read rank of each processor
    int rank;

    // Load each process rank
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)) { mpi_error(); }

    if (n_nodes == 1) {
        if (rank == ROOT_NODE) {

            std::cout << input_tree << std::endl;
        }
        if (MPI_Finalize()) { mpi_error(); }
        return EXIT_SUCCESS;
    }

    // Algorithms in PRL slides work with indexing from 1 -> follow convention
    unsigned edge_id = (unsigned) rank + 1;

    // Compute reverse edge id and weight
    unsigned reverse_edge = even(edge_id) ? edge_id - 1 : edge_id + 1;
    unsigned weight = (unsigned) even(rank);

    // Load adjacent representation, this is not part of pre-order algorithm, no need to optimize prize,
    // broadcasting from root node may be slower than calling with all processors
    load_adjacent_representation(target_nodes, next_edge, nodes_first_edge, tree_depth, n_nodes, n_edges);
    if (MPI_Barrier(MPI_COMM_WORLD)) { mpi_error(); }

//    // Debug print
//    if (rank == ROOT_NODE) {
//        prints_adjacent_representation(target_nodes, next_edge, nodes_first_edge);
//    }


    /* Euler tour
     * Get pointer to the next node in constant time by reading data from adjacent matrix
     * */

    unsigned e_next;

    if (next_edge[reverse_edge - 1] == NIL) {
        e_next = nodes_first_edge[target_nodes[rank] - 1];
    } else {
        e_next = next_edge[reverse_edge - 1];
    }

    // Parallel suffix sum presupposes last edge to contain cycle, process correction
    if (edge_id == LAST_EDGE) {
        e_next = edge_id;
    }

    /* Parallel suffix sum
     * No need to process first step by setting val of last edge to 0 and correction step,
     * because last edge is backward edge, so it's weight is already 0.
     * */
    std::vector<unsigned> succesors(n_edges);
    std::vector<unsigned> weights(n_edges);

    // Save original weight to filter out backward edges
    unsigned weight_original = weight;

    // Steps to be processed (log n complexity)
    unsigned pss_steps = ceil(log2((float) n_edges));

    if (MPI_Allgather(&e_next, 1, MPI_UNSIGNED, succesors.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD)) { mpi_error(); }
    if (MPI_Allgather(&weight, 1, MPI_UNSIGNED, weights.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD)) { mpi_error(); }

//    // Debug print
//    if (rank == ROOT_NODE) {
//        print_single_arr(succesors, "e_tour");
//    }

    // Core computation of suffix sum
    for (int i = 0; i < pss_steps; ++i) {
        weight += weights[e_next - 1];
        e_next = succesors[e_next - 1];
        if (MPI_Allgather(&e_next, 1, MPI_UNSIGNED, succesors.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD)) { mpi_error(); }
        if (MPI_Allgather(&weight, 1, MPI_UNSIGNED, weights.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD)) { mpi_error(); }
    }

    /* Correction step
     * */
    unsigned preorder = 0;
    if (weight_original) {
        preorder = n_nodes - weight;
    }

    /* Collect edge preorder ranking and print traversal to stdout
     * */
    std::vector<unsigned> preorders(n_edges);
    if (MPI_Gather(&preorder, 1, MPI_UNSIGNED, preorders.data(), 1, MPI_UNSIGNED, ROOT_NODE,
                   MPI_COMM_WORLD)) { mpi_error(); }

    if (rank == ROOT_NODE) {
//        // Debug print
//        print_single_arr(preorders, "preorders");
        std::string result(input_tree.length(), input_tree[0]);
        for (int i = 0; i < n_edges; i++) {
            unsigned edge = preorders[i];
            if (edge > 0) {
                unsigned input_position = target_nodes[i] - 1;
                result[edge] = input_tree[input_position];
            }
        }
        std::cout << result << std::endl;
    }

    if (MPI_Finalize()) { mpi_error(); }

    return EXIT_SUCCESS;
}