#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/property_map/property_map.hpp>

//===========================================================================
//// Base graph with custom info
//===========================================================================

struct VertexInfo {
  VertexInfo(const VertexInfo& vinfo) : qvi(vinfo.qvi) {}
  VertexInfo() : qvi(-1) {}
  int qvi;
};

//---------------------------------------------------------------------------

struct EdgeInfo {

  EdgeInfo(const EdgeInfo& einfo)
      : qei(einfo.qei),               //
        base_label(einfo.base_label), //
        curve(einfo.curve),           //
        length(einfo.length) {}

  EdgeInfo(int _qei, int _base_label, int _curve, double _length)
      : qei(_qei),               //
        base_label(_base_label), //
        curve(_curve),           //
        length(_length) {}

  EdgeInfo() : EdgeInfo(-1, -1, -1, 0.0) {}

  int qei;
  int base_label;
  int curve;
  double length; // for weights in dijkstra_shortest_paths via EdgeInfo::length
};

//---------------------------------------------------------------------------

using BaseGraph  = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeInfo>;
using BaseVertex = boost::graph_traits<BaseGraph>::vertex_descriptor;
using BaseSize   = boost::graph_traits<BaseGraph>::vertices_size_type;
using BaseEdge   = boost::graph_traits<BaseGraph>::edge_descriptor;
using QEdgeMap   = boost::property_map<BaseGraph, int EdgeInfo::*>::type;

//===========================================================================
//// Filtered graphs
//===========================================================================

struct Label_1 {
  Label_1() {}
  Label_1(int _sbl, QEdgeMap _ebl) : selected_base_label(_sbl), edge_base_labels(_ebl) {}
  QEdgeMap edge_base_labels;
  int selected_base_label;
  bool operator()(const BaseEdge& e) const { return edge_base_labels[e] == selected_base_label; }
};

//---------------------------------------------------------------------------

struct Label_0 {
  Label_0() {}
  Label_0(int _sbl, QEdgeMap _ebl) : selected_base_label(_sbl), edge_base_labels(_ebl) {}
  QEdgeMap edge_base_labels;
  int selected_base_label;
  bool operator()(const BaseEdge& e) const { return edge_base_labels[e] != selected_base_label; }
};

//---------------------------------------------------------------------------

using FilterGraph__Label_1 = boost::filtered_graph<BaseGraph, Label_1>; // predicate: edge has label
using FilterGraph__Label_0 = boost::filtered_graph<BaseGraph, Label_0>; // predicate: edge does not have label

//---------------------------------------------------------------------------
