#include <igsv/extraction/init_chains.h>
#include <igsv/common/defs.h>
#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/Graph.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QPort.h>
#include <igsv/extraction_entities/QVertex.h>

#include <set>

//===========================================================================
//// Depth-first search
//===========================================================================

class ChainVisitor : public boost::default_dfs_visitor {
public:
  IGSV::Chain& chain;
  ChainVisitor(IGSV::Chain& _chain) : chain(_chain) {}

  //// Visitor event points : https://www.boost.org/doc/libs/1_67_0/libs/graph/doc/depth_first_search.html
  // -- discover_vertex(u, g): invoked when a vertex is encountered for the first time
  void discover_vertex(BaseVertex v, const FilterGraph__Label_1& g) {
    const int deg = boost::degree(v, g);
    if (deg == 0) // ignore isolated vertices
      return;
    const int i = chain.qverts.size();
    chain.qverts.push_back(v);
    chain.qverts_degree.push_back(deg);
    chain.vindex_map.insert(std::make_pair(v, i));
    // chain.qverts2qedge.insert()
  }
  // -- finish_edge(e, g): invoked on the non-tree edges and on each tree edge after its target vertex is finished
  void finish_edge(BaseEdge e, const FilterGraph__Label_1& g) { //
    chain.qedges.insert(g[e].qei);
  }
};

//=============================================================================

namespace IGSV {

  void get_leaf_vertices(const FilterGraph__Label_1& _graph,   //
                         std::vector<BaseVertex>& _val1_verts, //
                         std::vector<BaseVertex>& _val2_verts) {

    _val1_verts.clear();
    _val2_verts.clear();

    FilterGraph__Label_1::vertex_iterator vit, vend;

    for (boost::tie(vit, vend) = boost::vertices(_graph); vit != vend; ++vit)
      switch (boost::degree(*vit, _graph)) {
      case 1:
        _val1_verts.push_back(*vit);
        break;

      case 2:
        _val2_verts.push_back(*vit);
        break;

      default:
        break;
      }

    std::sort(_val1_verts.begin(), _val1_verts.end());
    std::sort(_val2_verts.begin(), _val2_verts.end());

    _val1_verts.erase(std::unique(_val1_verts.begin(), _val1_verts.end()), _val1_verts.end());
    _val2_verts.erase(std::unique(_val2_verts.begin(), _val2_verts.end()), _val2_verts.end());
  }

  //===========================================================================

  bool get_two_most_distant_vertices(const FilterGraph__Label_1& _graph, //
                                     BaseVertex& _source, BaseVertex& _target) {

    std::vector<BaseVertex> leaves, _val2_verts;
    get_leaf_vertices(_graph, leaves, _val2_verts);
    if (leaves.size() < 2) {
      if (_val2_verts.empty())
        return false;

      // closed chain: source == target
      _source = (BaseVertex)_val2_verts.front();
      _target = _source;

      return true;
    }

    std::vector<double> distances(boost::num_vertices(_graph));

    //// distances between all pairs of leaf vertices
    Eigen::MatrixXd D;
    D.setZero(leaves.size(), leaves.size());

    int i = 0, j = 0;
    for (auto leaf : leaves) {

      boost::dijkstra_shortest_paths(_graph, leaf,
                                     boost::weight_map(boost::get(&EdgeInfo::length, _graph))
                                         .distance_map(boost::make_iterator_property_map(
                                             distances.begin(), boost::get(boost::vertex_index, _graph))));
      for (j = 0; j < leaves.size(); ++j)
        D(i, j) = distances[leaves[j]];

      ++i;
    }

    // source, target: which two leaf vertices are most distant? (longest curve)
    int rowId, colId;
    const double max = D.maxCoeff(&rowId, &colId);

    _source = (BaseVertex)leaves[rowId];
    _target = (BaseVertex)leaves[colId];

    return true;
  }

  //===========================================================================

  void get_adjacent_edges(FilterGraph__Label_1& g_label, //
                          FilterGraph__Label_0& g_compl, //
                          std::set<int>& _edges) {

    // loop over edges of the complement
    FilterGraph__Label_0::edge_iterator eit, eend;
    for (boost::tie(eit, eend) = boost::edges(g_compl); eit != eend; ++eit) {
      BaseVertex p = boost::source(*eit, g_compl); // source vertex
      BaseVertex q = boost::target(*eit, g_compl); // target vertex
      const int dp = boost::degree(p, g_label);    // degree of p in g_label
      const int dq = boost::degree(q, g_label);    // degree of q in g_label
      if (dp > 0 || dq > 0) {                      // one of the vertices is in g_label
        _edges.insert(g_compl[*eit].qei);
      }
    }
  }

  //===========================================================================

  bool init_chains(const std::vector<QVertex>& _QV, //
                   std::vector<QEdge>& _QE,         //
                   std::vector<Chain>& _CH) {

    _CH.clear();

    const int n_qverts = _QV.size();
    const int n_qedges = _QE.size();

    int chid = 0;

    // ===========================================================================
    //// Construct graph
    // ===========================================================================
    BaseGraph base_graph(n_qverts);
    std::set<int> base_labels;

    // add vertices
    for (int qvi = 0; qvi < n_qverts; ++qvi)
      base_graph[qvi].qvi = qvi;

    // add edges
    for (int qei = 0; qei < n_qedges; ++qei)
      if (_QE[qei].base_label > -1) {
        boost::add_edge(_QE[qei].qvi0, _QE[qei].qvi1, EdgeInfo(qei, _QE[qei].base_label, -1, _QE[qei].length),
                        base_graph);
        base_labels.insert(_QE[qei].base_label);
      }

    // ===========================================================================
    //// Construct chains
    // ===========================================================================

    const QEdgeMap edge_base_labels = boost::get(&EdgeInfo::base_label, base_graph);

    int cid = 0; // chain index

    for (auto base_label : base_labels) {

      FilterGraph__Label_1 g_label(base_graph, Label_1(base_label, edge_base_labels));
      FilterGraph__Label_0 g_compl(base_graph, Label_0(base_label, edge_base_labels));

      // get vertex valences for the filtered graph
      std::vector<int> valences(n_qverts, 0);
      int n_fgraph_edges = 0;
      FilterGraph__Label_1::edge_iterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(g_label); eit != eend; ++eit) {
        valences[boost::source(*eit, g_label)]++;
        valences[boost::target(*eit, g_label)]++;
        n_fgraph_edges++;
      }

      BaseVertex source, target;
      if (!get_two_most_distant_vertices(g_label, source, target)) {
        DEBUG_PRINT_WARNING("chain #%d: get_two_most_distant_vertices returned false", cid);
        continue;
      }

      Chain chain;
      chain.qvi0 = (int)source;
      chain.qvi1 = (int)target;
      ChainVisitor chain_visitor(chain);
      boost::depth_first_search(g_label, boost::visitor(chain_visitor).root_vertex(chain.qvi0));

      get_adjacent_edges(g_label, g_compl, chain.adjacent_qedges);

      for (auto qei : chain.qedges)
        _QE[qei].chain = chid;

      _CH.push_back(chain);
      chid++;

    } // end loop over chains

    DEBUG_PRINT_INFO("%zd chains", _CH.size());

    return true;
  }

} // namespace IGSV
