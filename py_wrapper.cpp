#include <pybind11/pybind11.h>
#include <iostream>

// declaration of the graph type lives here
#include "graphs_common.h"

#include "compute.h"

namespace py = pybind11;

void process(py::object nx_graph, long long time_limit) {
    // import networkx in C++ to access helper functions
    py::module_ nx = py::module_::import("networkx");

    // convert nodes to 0-indexed integers safely (handles string nodes, etc.)
    py::object clean_graph = nx.attr("convert_node_labels_to_integers")(nx_graph);

    // allocate the graph_data struct on the heap
    graph g = new graph_data();

    // basic properties
    g->n = clean_graph.attr("number_of_nodes")().cast<refer>();
    g->m = clean_graph.attr("number_of_edges")().cast<unsigned long>();
    g->density = nx.attr("density")(clean_graph).cast<double>();

    // safety check (until we get rid of MAX_VERTICES)
    if (g->n > MAX_VERTICES) {
        delete g;
        throw std::runtime_error("Graph exceeds MAX_VERTICES");
    }

    // extract adjacency list (dictionary of lists in Python)
    py::dict adj = nx.attr("to_dict_of_lists")(clean_graph);

    for (int i = 0; i < g->n; ++i) {
        py::list neighbors = adj[py::cast(i)].cast<py::list>();
        g->V[i].edgecount = neighbors.size();

        // memory allocation for siblings
        g->V[i].sibl = new refer[g->V[i].edgecount];

        int j = 0;
        for (auto neighbor : neighbors) {
            g->V[i].sibl[j++] = neighbor.cast<refer>();
        }
    }

//    std::cout << "Nodes: " << g->n << ", Edges: " << g->m << ", Density: " << g->density << "\n";
//    for (int i = 0; i < g->n; ++i) {
//        std::cout << "Node " << i << " -> ";
//        for (int j = 0; j < g->V[i].edgecount; ++j) {
//            std::cout << g->V[i].sibl[j] << " ";
//        }
//        std::cout << "\n";
//    }

    refer *coloring = new refer[g->n];

    // call the solver
    compute(g, coloring, time_limit);

    delete[](coloring);

    // cleanup
    for (int i = 0; i < g->n; ++i) {
        delete[] g->V[i].sibl;
    }
    delete g;
}

PYBIND11_MODULE(coloring_utility, m) {
    m.doc() = "Python wrapper for ColoringUtility";
    // export of our wrapper function
    m.def(
        "process",
        &process,
        "Accepts a networkx graph and computes a coloring using the ColoringUtility solver",
        py::arg("nx_graph"),
        py::arg("time_limit") = 60
    );
}
