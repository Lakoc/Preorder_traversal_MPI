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

void prints_adjacent_representation(std::vector<unsigned> weights, std::vector<unsigned> target_nodes,
                                    std::vector<int> next_edge, std::vector<unsigned> nodes_first_edge) {
    debug_print("weights: ");
    for (unsigned i: weights) {
        debug_print(i << " ");
    }
    debug_print("\ntarget_nodes: ");
    for (unsigned i: target_nodes) {
        debug_print(i << " ");
    }
    debug_print("\nnext_edge: ");
    for (int i: next_edge) {
        debug_print(i << " ");
    }
    debug_print("\nnodes_first_edge: ");
    for (int i: nodes_first_edge) {
        debug_print(i << " ");
    }
    debug_print("\n");
}


void load_adjacent_representation(std::vector<unsigned>& weights, std::vector<unsigned>& target_nodes,
                                  std::vector<int>& next_edge, std::vector<unsigned>& nodes_first_edge,
                                  unsigned tree_depth, unsigned n_nodes, unsigned n_edges) {
    unsigned first_forward = ROOT_NODE + 1;
    unsigned backward_start = ROOT_NODE * N_CHILDREN;
    unsigned even_offset = EVEN_OFFSET_START;
    unsigned offsets_counted = 0;

    for (unsigned i = 0; i < n_edges; ++i) {
        if (even(i)) {
            weights.push_back(1);
            target_nodes.push_back(first_forward);
            first_forward += 1;
            if (even(i / 2)) {
                next_edge.push_back(N_CHILDREN);
            } else {
                next_edge.push_back(NILL);
            }
        } else {
            weights.push_back(0);
            target_nodes.push_back(backward_start / BACKWARD_NODES);
            backward_start += 1;

            if (offsets_counted < tree_depth) {
                next_edge.push_back(even_offset);
                even_offset += N_CHILDREN;
                offsets_counted += 1;
            } else {
                next_edge.push_back(NILL);
            }
        }
    }

    nodes_first_edge.push_back(ROOT_NODE);
    for (unsigned i = 1; i < n_nodes; ++i) {
        nodes_first_edge.push_back(i * N_CHILDREN);
    }
}

int main(int argc, char **argv) {

    int n_processes;
    int group, edge_id, reverse_edge;
    std::string input_tree = argv[1]; // Name of the current exec program
    unsigned n_nodes = input_tree.length();
    unsigned n_edges = (2 * n_nodes) - 2;
    unsigned tree_depth = log2base(n_nodes);

    std::vector<unsigned> weights;
    std::vector<unsigned> target_nodes;
    std::vector<int> next_edge;

    std::vector<unsigned> nodes_first_edge;


    if (MPI_Init(&argc, &argv)) { mpi_error(); }

    // Load process counts and each process rank
    if (MPI_Comm_size(MPI_COMM_WORLD, &n_processes)) { mpi_error(); }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &edge_id)) { mpi_error(); }

    edge_id += 1;
    reverse_edge = even(edge_id) ? edge_id - 1 : edge_id + 1;

    load_adjacent_representation(weights, target_nodes, next_edge, nodes_first_edge, tree_depth, n_nodes, n_edges);

    if (edge_id == ROOT_NODE) {
        if (DEBUG) {
            prints_adjacent_representation(weights, target_nodes, next_edge, nodes_first_edge);
        }

    }

    if (MPI_Barrier(MPI_COMM_WORLD)) { mpi_error(); }

    if (MPI_Finalize()) { mpi_error(); }

    return EXIT_SUCCESS;
}