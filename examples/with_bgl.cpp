// Examples of using Boost Graph Library with gemmi::ChemComp:
//  checking graph and subgraph isomorphism,
//  finding maximum common induced subgraph.
//
// NOTE: documentation (mol.rst) includes parts of this file.

#include <cstring>                   // for strcmp
#include <iostream>                  // for std::cout, std::cerr
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>
#include <gemmi/cif.hpp>             // for cif::read_file
#include <gemmi/chemcomp.hpp>        // for ChemComp, make_chemcomp_from_block

struct AtomVertex {
  gemmi::Element el = gemmi::El::X;
  std::string name;
};

struct BondEdge {
  gemmi::BondType type;
};

using Graph = boost::adjacency_list<boost::vecS, boost::vecS,
                                    boost::undirectedS,
                                    AtomVertex, BondEdge>;

Graph make_graph(const gemmi::ChemComp& cc) {
  Graph g(cc.atoms.size());
  for (size_t i = 0; i != cc.atoms.size(); ++i) {
    g[i].el = cc.atoms[i].el;
    g[i].name = cc.atoms[i].id;
  }
  for (const gemmi::Restraints::Bond& bond : cc.rt.bonds) {
    int n1 = cc.get_atom_index(bond.id1.atom);
    int n2 = cc.get_atom_index(bond.id2.atom);
    boost::add_edge(n1, n2, BondEdge{bond.type}, g);
  }
  return g;
}

gemmi::ChemComp make_chemcomp(const char* path) {
  gemmi::cif::Document doc = gemmi::cif::read_file(path);
  // assuming the component description is in the last block of the file
  return gemmi::make_chemcomp_from_block(doc.blocks.back());
}

using GraphTraits = boost::graph_traits<Graph>;
using Vertex = GraphTraits::vertex_descriptor;
using Edge = GraphTraits::edge_descriptor;


// Example 1 - count automorphisms (a toy example).
struct CountingCallback {
  int n = 0;
  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) {
    ++n;
    return true;
  }
};

void count_automorphisms_of_SO3() {
  Graph g = make_graph(make_chemcomp("tests/SO3.cif"));
  CountingCallback counting_callback;
  boost::vf2_graph_iso(
          g, g, std::ref(counting_callback),
          // default values
          boost::get(boost::vertex_index, g),
          boost::get(boost::vertex_index, g),
          boost::vertex_order_by_mult(g),
          // ignoring bond types - all edges are considered equivalent
          boost::always_equivalent(),
          // atoms of the same type (element) are considered equivalent
          [&](Vertex a, Vertex b) { return g[a].el == g[b].el; });
  std::cout << counting_callback.n << std::endl;
  // prints 6 (3! permutations of the three oxygens in SO3)
}


// Example 2 - check graph isomorphism.
//
// Example output:
//
//  $ ./with_bgl ccd/M10.cif monomers/m/M10.cif 
//  isomorphic!
//    O4 -> O9
//    O9 -> O4
//
struct PrintMappingCallback {
  Graph g1, g2;
  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 map1, CorrespondenceMap2To1) {
    std::cout << "isomorphic!\n";
    for (auto vp = boost::vertices(g1); vp.first != vp.second; ++vp.first) {
      Vertex v1 = *vp.first;
      Vertex v2 = boost::get(map1, v1);
      if (v2 != GraphTraits::null_vertex()) {
        const std::string& name1 = g1[v1].name;
        const std::string& name2 = g2[v2].name;
        if (name1 != name2)
          std::cout << "  " << name1 << " -> " << name2 << std::endl;
      }
    }
    return false;
  }
};

void check_graph_isomorphism(const char* cif1, const char* cif2) {
  Graph graph1 = make_graph(make_chemcomp(cif1));
  Graph graph2 = make_graph(make_chemcomp(cif2));
  bool found = vf2_graph_iso(
      graph1, graph2,
      PrintMappingCallback{graph1, graph2},
      get(boost::vertex_index, graph1), get(boost::vertex_index, graph2),
      boost::vertex_order_by_mult(graph1),
      [&](Edge a, Edge b) { return graph1[a].type == graph2[b].type; },
      [&](Vertex a, Vertex b) { return graph1[a].el == graph2[b].el; });
  if (!found)
    std::cout << "not isomorphic" << std::endl;
}


// Example 3 - check substructure isomorphism (if any),
// ignoring hydrogens and bond types.
void check_subgraph_isomorphism(const char* cif1, const char* cif2) {
  gemmi::ChemComp cc1 = make_chemcomp(cif1).remove_hydrogens();
  Graph graph1 = make_graph(cc1);
  gemmi::ChemComp cc2 = make_chemcomp(cif2).remove_hydrogens();
  Graph graph2 = make_graph(cc2);
  bool found = vf2_subgraph_iso(
      graph1, graph2,
      PrintMappingCallback{graph1, graph2},
      get(boost::vertex_index, graph1), get(boost::vertex_index, graph2),
      boost::vertex_order_by_mult(graph1),
      boost::always_equivalent(), // edge_comp
      [&](Vertex a, Vertex b) { return graph1[a].el == graph2[b].el; });
  if (!found)
    std::cout << "not isomorphic" << std::endl;
}


// Example 4 - find the largest common subgraph using McGregor's algorithm.
// We ignore hydrogens here.
// Example output:
//   $ time with_bgl -c monomers/a/AUD.cif monomers/l/LSA.cif
//   Searching largest subgraphs of AUD and LSA (10 and 12 vertices)...
//   Maximum connected common subgraph has 7 vertices:
//     SAA -> S10
//     OAG -> O12
//     OAH -> O11
//     NAF -> N9
//     CAE -> C7
//     OAI -> O8
//     CAD -> C1
//   Maximum connected common subgraph has 7 vertices:
//     SAA -> S10
//     OAG -> O11
//     OAH -> O12
//     NAF -> N9
//     CAE -> C7
//     OAI -> O8
//     CAD -> C1
//   
//   real	0m0.012s
//   user	0m0.008s
//   sys	0m0.004s

struct PrintCommonSubgraphCallback {
  Graph g1, g2;
  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 map1, CorrespondenceMap2To1,
                  typename GraphTraits::vertices_size_type subgraph_size) {
    std::cout << "Maximum connected common subgraph has " << subgraph_size
              << " vertices:\n";
    for (auto vp = boost::vertices(g1); vp.first != vp.second; ++vp.first) {
      Vertex v1 = *vp.first;
      Vertex v2 = boost::get(map1, v1);
      if (v2 != GraphTraits::null_vertex())
        std::cout << "  " << g1[v1].name << " -> " << g2[v2].name << std::endl;
    }
    return true;
  }
};

void find_common_subgraph(const char* cif1, const char* cif2) {
  gemmi::ChemComp cc1 = make_chemcomp(cif1).remove_hydrogens();
  Graph graph1 = make_graph(cc1);
  gemmi::ChemComp cc2 = make_chemcomp(cif2).remove_hydrogens();
  Graph graph2 = make_graph(cc2);
  std::cout << "Searching largest subgraphs of " << cc1.name << " and "
      << cc2.name << " (" << cc1.atoms.size() << " and " << cc2.atoms.size()
      << " vertices)..." << std::endl;
  mcgregor_common_subgraphs_maximum_unique(
      graph1, graph2,
      get(boost::vertex_index, graph1), get(boost::vertex_index, graph2),
      [&](Edge a, Edge b) { return graph1[a].type == graph2[b].type; },
      [&](Vertex a, Vertex b) { return graph1[a].el == graph2[b].el; },
      /*only_connected_subgraphs=*/ true,
      PrintCommonSubgraphCallback{graph1, graph2});
}


// a minimal program that can call the functions above
int main(int argc, char* argv[]) {
  if (argc == 1)
    count_automorphisms_of_SO3(); // toy example that outputs "6".
  else if (argc == 3)
    check_graph_isomorphism(argv[1], argv[2]);
  else if (argc == 4 && std::strcmp(argv[1], "-s") == 0)
    check_subgraph_isomorphism(argv[2], argv[3]);
  else if (argc == 4 && std::strcmp(argv[1], "-c") == 0)
    find_common_subgraph(argv[2], argv[3]);
  else
    std::cerr << "Usage: " << argv[0] << " [-s|-c] ligand1.cif ligand2.cif\n";
  return 0;
}
