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
 * Create adjacent matrix in constant time parallel for each processor
 *
 *  target_node: destination of current edge
 *
 *  next_edge: `pointer` to the next edge, contains NIL if next edge does not exist
 *
 *  tree_depth: depth of tree - 1
 *
 *  n_edges: number of edges in tree
 */
void
load_adjacent_representation(int rank, unsigned *target_node, int *next_edge,
                             unsigned tree_depth, unsigned n_edges) {
    unsigned edge_id = (unsigned) rank + 1;
    if (even(rank)) {
        *target_node = ROOT_ID + 1 + (rank / 2);
        if (even(rank / 2) && edge_id + N_CHILDREN <= n_edges) {
            *next_edge = edge_id + N_CHILDREN;
        } else {
            *next_edge = NIL;
        }

    } else {
        *target_node = ((rank / 2) + (N_CHILDREN)) / BACKWARD_NODES;
        if (edge_id + (rank / 2 * N_CHILDREN) + EVEN_OFFSET_START <= n_edges) {
            if (tree_depth <= 1) {
                *next_edge = (NIL);
            } else {
                *next_edge = edge_id + (rank / 2 * N_CHILDREN) + EVEN_OFFSET_START;
            }
        } else {
            *next_edge = NIL;
        }
    }
}

/*
 * Function:  get_reversed_next_edge
 * --------------------
 * Request data from reversed processor
 *
 *  rank: curr processor rank
 *
 *  reversed_rank: reversed edge processor rank
 *
 *  next_edge: processor data to be sent
 *
 *  returns: reversed processor next edge
 */
int get_reversed_next_edge(unsigned rank, unsigned reversed_rank, int next_edge) {
    MPI_Request req;
    int reversed_next_edge;

    if (MPI_Isend(&next_edge, 1, MPI_INT, reversed_rank, DEFAULT_TAG, MPI_COMM_WORLD, &req)) { mpi_error(); }
    if (MPI_Recv(&reversed_next_edge, 1, MPI_INT, MPI_ANY_SOURCE, DEFAULT_TAG, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE)) { mpi_error(); }

    if (MPI_Request_free(&req)) { mpi_error(); }
    return reversed_next_edge;
}

int main(int argc, char **argv) {
    if (MPI_Init(&argc, &argv)) { mpi_error(); }
    // Load input tree
    std::string input_tree = argv[1]; // Name of the current exec program

    // Calculate number of edges, nodes and tree depth
    unsigned n_nodes = input_tree.length();
    unsigned n_edges = (2 * n_nodes) - 2;
    unsigned tree_depth = ceil(log2(n_nodes));

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

    // Load adjacent representation, compute for each edge target node and next edge in constant time
    unsigned target_node;
    int next_edge;
    load_adjacent_representation(rank, &target_node, &next_edge, tree_depth, n_edges);

    /* Euler tour
     * Get pointer to the next node in constant time by reading data from adjacent matrix
     * */

    unsigned e_next;

    int reversed_next_rank = get_reversed_next_edge(rank, reverse_edge - 1, next_edge);

    if (reversed_next_rank == NIL) {
        e_next = target_node == 1 ? target_node : (target_node - 1) * 2;
    } else {
        e_next = reversed_next_rank;
    }

    /* Parallel suffix sum
     * No need to process first step by setting val of last edge to 0 and correction step,
     * because last edge is backward edge, so it's weight is already 0.
     * */

    // Parallel suffix sum presupposes last edge to contain cycle, process correction
    if (edge_id == LAST_EDGE) {
        e_next = edge_id;
    }

    std::vector<unsigned> succesors(n_edges);
    std::vector<unsigned> weights(n_edges);

    // Save original weight to filter out backward edges
    unsigned weight_original = weight;

    // Steps to be processed (log n complexity)
    unsigned pss_steps = ceil(log2((float) n_edges));

    // Synchronize information between nodes
    if (MPI_Allgather(&e_next, 1, MPI_UNSIGNED, succesors.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD)) { mpi_error(); }
    if (MPI_Allgather(&weight, 1, MPI_UNSIGNED, weights.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD)) { mpi_error(); }

//    // Debug print
//    if (rank == ROOT_NODE) {
//        print_single_arr(succesors, "e_tour");
//    }

    // Core computation of suffix sum
    for (int i = 0; i < pss_steps; i++) {
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
    std::vector<unsigned> target_nodes(n_edges);

    if (MPI_Gather(&preorder, 1, MPI_UNSIGNED, preorders.data(), 1, MPI_UNSIGNED, ROOT_NODE,
                   MPI_COMM_WORLD)) { mpi_error(); }
    if (MPI_Gather(&target_node, 1, MPI_UNSIGNED, target_nodes.data(), 1, MPI_UNSIGNED, ROOT_NODE,
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