#include <igsv/extraction/connected_components.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <map>

namespace IGSV {
  // based on https://stackoverflow.com/questions/46659904/connected-components-in-an-undirected-graph-in-c
  void connected_components(const Eigen::MatrixXi& E, std::vector<std::vector<int>>& conncomp, int n) {
    assert(E.cols() == 2);
    using Graph   = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
    using Mapping = std::map<Graph::vertex_descriptor, int>;

    // construct the graph
    Graph graph;
    for (int v = 0; v < n; v++) {
      boost::add_edge(v, v, graph);
    }
    for (int e = 0; e < E.rows(); e++) {
      boost::add_edge(E(e, 0), E(e, 1), graph);
    }

    // compute the connected components
    Mapping mappings;
    const int nc = boost::connected_components(graph, boost::make_assoc_property_map(mappings));

    // store
    conncomp.resize(nc);
    for (auto& mapping : mappings) {
      conncomp[mapping.second].push_back(mapping.first);
    }
  }

  void connected_components(const std::vector<std::pair<int, int>>& edges, //
                            std::vector<std::vector<int>>& conncomp,       //
                            int n) {
    using Graph   = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
    using Mapping = std::map<Graph::vertex_descriptor, int>;

    // construct the graph
    Graph graph;
    for (int v = 0; v < n; v++) {
      boost::add_edge(v, v, graph);
    }
    for (const auto e : edges) {
      boost::add_edge(e.first, e.second, graph);
    }

    // compute the connected components
    Mapping mappings;
    const int nc = boost::connected_components(graph, boost::make_assoc_property_map(mappings));

    // store
    conncomp.resize(nc);
    for (auto& mapping : mappings) {
      conncomp[mapping.second].push_back(mapping.first);
    }
  }

}; // namespace IGSV
