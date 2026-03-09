// Copyright 2026 Global Phasing Ltd.
// Internal carborane processing for AceDRG-style pipeline.

#include "gemmi/ace_carborane.hpp"
#include "gemmi/ace_graph.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/fileutil.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace gemmi {

namespace {

template<typename Pred>
std::set<size_t> collect_carborane_atoms(const ChemComp& cc,
                                         const AceBondAdjacency& adj,
                                         Pred should_expand) {
  std::vector<size_t> seeds;
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    if (cc.atoms[i].el == El::H)
      continue;
    int b_count = 0;
    for (const auto& nb : adj[i])
      if (cc.atoms[nb.idx].el == El::B)
        ++b_count;
    if (b_count >= 4)
      seeds.push_back(i);
  }
  std::set<size_t> out;
  std::vector<size_t> stack = seeds;
  while (!stack.empty()) {
    size_t idx = stack.back();
    stack.pop_back();
    if (!out.insert(idx).second)
      continue;
    for (const auto& nb : adj[idx]) {
      Element el = cc.atoms[nb.idx].el;
      if (el != El::H && should_expand(el, nb.idx))
        stack.push_back(nb.idx);
    }
  }
  return out;
}

std::set<size_t> collect_carborane_cluster_atoms(const ChemComp& cc,
                                                  const AceBondAdjacency& adj) {
  return collect_carborane_atoms(cc, adj, [](Element el, size_t) {
    return el == El::B || el == El::C || el.is_metal();
  });
}

std::set<size_t> collect_carborone_match_atoms(const ChemComp& cc,
                                                const AceBondAdjacency& adj) {
  std::vector<int> b_neighbors(cc.atoms.size(), 0);
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    for (const auto& nb : adj[i])
      if (cc.atoms[nb.idx].el == El::B)
        ++b_neighbors[i];
  return collect_carborane_atoms(cc, adj, [&](Element el, size_t idx) {
    return el == El::B || el.is_metal() || b_neighbors[idx] >= 2;
  });
}

struct CarboroneGraph {
  std::vector<std::string> atom_ids;
  std::vector<std::string> elem_syms;
  std::vector<Position> xyz;
  std::vector<std::vector<int>> edge_w;
  std::vector<std::vector<int>> neighbors;
};

struct CarboroneDb {
  bool list_loaded = false;
  std::map<std::string, std::vector<std::string>> by_pair;
  std::map<std::string, std::vector<std::string>> by_class0;
  std::vector<std::string> all_names;
  std::map<std::string, CarboroneGraph> templates;
};

std::string upper_copy(std::string s) {
  for (char& c : s)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

std::string join_path(const std::string& dir, const std::string& name) {
  if (dir.empty())
    return name;
  if (dir.back() == '/')
    return dir + name;
  return dir + "/" + name;
}

int carborone_weight_from_symbol(const std::string& sym_u) {
  if (sym_u == "B") return 2;
  if (sym_u == "C") return 3;
  if (sym_u == "N") return 4;
  if (sym_u == "O") return 5;
  if (sym_u == "S") return 6;
  if (sym_u == "SE") return 7;
  if (sym_u == "P") return 8;
  if (sym_u == "AS" || sym_u == "SI" || sym_u == "GA" ||
      sym_u == "GI" || sym_u == "IN")
    return 9;
  return 1;
}

void build_carborone_graph(const ChemComp& cc,
                           const std::vector<size_t>& atom_indices,
                           CarboroneGraph& out,
                           std::vector<int>* global_to_local) {
  out = CarboroneGraph();
  std::vector<int> idx_to_local(cc.atoms.size(), -1);
  for (size_t idx : atom_indices) {
    if (idx >= cc.atoms.size() || cc.atoms[idx].el == El::H)
      continue;
    int local = static_cast<int>(out.atom_ids.size());
    idx_to_local[idx] = local;
    out.atom_ids.push_back(cc.atoms[idx].id);
    out.elem_syms.push_back(upper_copy(cc.atoms[idx].el.uname()));
    out.xyz.push_back(cc.atoms[idx].xyz);
  }
  size_t n = out.atom_ids.size();
  out.edge_w.assign(n, std::vector<int>(n, 0));
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom);
    int idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0) continue;
    int l1 = idx_to_local[(size_t) idx1];
    int l2 = idx_to_local[(size_t) idx2];
    if (l1 < 0 || l2 < 0 || l1 == l2) continue;
    int w = carborone_weight_from_symbol(out.elem_syms[(size_t) l1]) +
            carborone_weight_from_symbol(out.elem_syms[(size_t) l2]);
    out.edge_w[(size_t) l1][(size_t) l2] = w;
    out.edge_w[(size_t) l2][(size_t) l1] = w;
  }
  out.neighbors.assign(n, std::vector<int>());
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      if (out.edge_w[i][j] != 0)
        out.neighbors[i].push_back((int) j);
  if (global_to_local)
    *global_to_local = std::move(idx_to_local);
}

std::vector<std::string> ordered_carborone_elements(const std::map<std::string, int>& counts) {
  std::vector<std::string> out;
  if (counts.count("B")) out.push_back("B");
  if (counts.count("C")) out.push_back("C");
  for (const auto& kv : counts)
    if (kv.first != "B" && kv.first != "C")
      out.push_back(kv.first);
  return out;
}

void carborone_class_ids(const CarboroneGraph& graph,
                         std::string& class0,
                         std::string& class1) {
  std::map<std::string, int> elem_counts;
  std::map<std::string, std::map<int, int>> elem_conn_counts;
  for (size_t i = 0; i < graph.elem_syms.size(); ++i) {
    const std::string& sym = graph.elem_syms[i];
    elem_counts[sym] += 1;
    elem_conn_counts[sym][(int) graph.neighbors[i].size()] += 1;
  }
  std::vector<std::string> ordered = ordered_carborone_elements(elem_counts);
  class0.clear();
  class1.clear();
  for (const std::string& sym : ordered)
    class0 += cat("[", sym, elem_counts[sym], "]");
  for (const std::string& sym : ordered)
    for (const auto& kv : elem_conn_counts[sym])
      class1 += cat("[", sym, kv.first, "_", kv.second, "]");
}

CarboroneDb& get_carborone_db(const std::string& tables_dir) {
  static std::map<std::string, CarboroneDb> cache;
  CarboroneDb& db = cache[tables_dir];
  if (db.list_loaded) return db;
  db.list_loaded = true;
  std::string path = join_path(tables_dir, "CarboroneSamples.list");
  fileptr_t f = file_open_or_null(path.c_str(), "r");
  if (!f) return db;
  std::set<std::string> seen_names;
  char buf[256];
  while (std::fgets(buf, sizeof(buf), f.get())) {
    char c0[64], c1[64], nm[64];
    if (std::sscanf(buf, "%63s %63s %63s", c0, c1, nm) != 3) continue;
    std::string class0 = upper_copy(c0);
    std::string class1 = upper_copy(c1);
    std::string name(nm);
    db.by_pair[cat(class0, '\t', class1)].push_back(name);
    db.by_class0[class0].push_back(name);
    if (seen_names.insert(name).second)
      db.all_names.push_back(name);
  }
  return db;
}

const CarboroneGraph* get_carborone_template(CarboroneDb& db,
                                             const std::string& tables_dir,
                                             const std::string& name) {
  auto found = db.templates.find(name);
  if (found != db.templates.end())
    return found->second.atom_ids.empty() ? nullptr : &found->second;
  CarboroneGraph graph;
  std::string path = join_path(join_path(tables_dir, "CarboroneSamples"), name + ".cif");
  try {
    cif::Document doc = cif::read_file(path);
    if (!doc.blocks.empty()) {
      ChemComp tcc = make_chemcomp_from_block(doc.blocks[0]);
      std::vector<size_t> idx;
      for (size_t i = 0; i < tcc.atoms.size(); ++i)
        if (tcc.atoms[i].el != El::H)
          idx.push_back(i);
      build_carborone_graph(tcc, idx, graph, nullptr);
    }
  } catch (...) {}
  db.templates[name] = graph;
  auto it = db.templates.find(name);
  return it->second.atom_ids.empty() ? nullptr : &it->second;
}

std::vector<std::string> collect_carborone_candidates(const CarboroneDb& db,
                                                      const std::string& class0,
                                                      const std::string& class1) {
  std::vector<std::string> out;
  std::set<std::string> seen;
  auto append_unique = [&](const std::vector<std::string>& src) {
    for (const std::string& s : src)
      if (seen.insert(s).second)
        out.push_back(s);
  };
  auto it_pair = db.by_pair.find(class0 + "\t" + class1);
  if (it_pair != db.by_pair.end()) append_unique(it_pair->second);
  if (out.empty()) {
    auto it0 = db.by_class0.find(class0);
    if (it0 != db.by_class0.end()) append_unique(it0->second);
  }
  if (out.empty()) append_unique(db.all_names);
  return out;
}

bool backtrack_carborone_match(const CarboroneGraph& target,
                               const CarboroneGraph& templ,
                               const std::vector<std::vector<int>>& candidates,
                               const std::vector<int>& order,
                               int depth,
                               std::vector<int>& t2s,
                               std::vector<int>& s2t) {
  if (depth == (int) order.size()) return true;
  int t = order[(size_t) depth];
  for (int s : candidates[(size_t) t]) {
    if (s2t[(size_t) s] != -1) continue;
    bool ok = true;
    for (size_t t2 = 0; t2 < t2s.size(); ++t2) {
      int s2 = t2s[t2];
      if (s2 == -1) continue;
      if (target.edge_w[(size_t) t][t2] != templ.edge_w[(size_t) s][(size_t) s2]) {
        ok = false; break;
      }
    }
    if (!ok) continue;
    for (size_t u = 0; u < t2s.size() && ok; ++u) {
      if (t2s[u] != -1 || u == (size_t) t) continue;
      int w = target.edge_w[(size_t) t][u];
      if (w == 0) continue;
      bool found = false;
      for (int su : candidates[u]) {
        if (s2t[(size_t) su] != -1) continue;
        if (templ.edge_w[(size_t) s][(size_t) su] == w) {
          found = true; break;
        }
      }
      if (!found) ok = false;
    }
    if (!ok) continue;
    t2s[(size_t) t] = s; s2t[(size_t) s] = t;
    if (backtrack_carborone_match(target, templ, candidates, order, depth + 1, t2s, s2t))
      return true;
    t2s[(size_t) t] = -1; s2t[(size_t) s] = -1;
  }
  return false;
}

bool match_carborone_graphs(const CarboroneGraph& target,
                            const CarboroneGraph& templ,
                            std::vector<int>& t2s_out) {
  const size_t n = target.elem_syms.size();
  if (n == 0 || templ.elem_syms.size() != n) return false;
  std::vector<std::vector<int>> target_inc_w(n), templ_inc_w(n);
  for (size_t i = 0; i < n; ++i) {
    for (int nb : target.neighbors[i])
      target_inc_w[i].push_back(target.edge_w[i][(size_t) nb]);
    std::sort(target_inc_w[i].begin(), target_inc_w[i].end());
    for (int nb : templ.neighbors[i])
      templ_inc_w[i].push_back(templ.edge_w[i][(size_t) nb]);
    std::sort(templ_inc_w[i].begin(), templ_inc_w[i].end());
  }
  std::vector<std::vector<int>> candidates(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (target.elem_syms[i] != templ.elem_syms[j]) continue;
      if (target.neighbors[i].size() != templ.neighbors[j].size()) continue;
      if (target_inc_w[i] != templ_inc_w[j]) continue;
      candidates[i].push_back((int) j);
    }
    if (candidates[i].empty()) return false;
  }
  std::vector<int> order(n);
  for (size_t i = 0; i < n; ++i) order[i] = (int) i;
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    size_t ca = candidates[(size_t) a].size(), cb = candidates[(size_t) b].size();
    return ca != cb ? ca < cb : target.neighbors[(size_t) a].size() > target.neighbors[(size_t) b].size();
  });
  std::vector<int> t2s(n, -1), s2t(n, -1);
  if (!backtrack_carborone_match(target, templ, candidates, order, 0, t2s, s2t))
    return false;
  t2s_out = std::move(t2s);
  return true;
}

bool apply_carborone_template_bonds(ChemComp& cc,
                                    const std::set<size_t>& cb_atoms,
                                    const std::string& tables_dir) {
  if (tables_dir.empty() || cb_atoms.empty()) return false;
  std::vector<size_t> target_idx;
  for (size_t idx : cb_atoms)
    if (idx < cc.atoms.size() && cc.atoms[idx].el != El::H)
      target_idx.push_back(idx);
  if (target_idx.size() < 3) return false;
  CarboroneGraph target;
  build_carborone_graph(cc, target_idx, target, nullptr);
  if (target.elem_syms.empty()) return false;
  std::string class0, class1;
  carborone_class_ids(target, class0, class1);
  CarboroneDb& db = get_carborone_db(tables_dir);
  std::vector<std::string> candidates = collect_carborone_candidates(db, class0, class1);
  for (const std::string& name : candidates) {
    const CarboroneGraph* templ = get_carborone_template(db, tables_dir, name);
    std::vector<int> dummy;
    if (templ && match_carborone_graphs(target, *templ, dummy))
      return true;
  }
  return false;
}

bool is_carborane_h_center(const ChemComp& cc, const AceBondAdjacency& adj, size_t idx) {
  const Element el = cc.atoms[idx].el;
  if (el != El::B && el != El::C) return false;
  bool has_h = false;
  int non_h = 0;
  for (const auto& nb : adj[idx]) {
    const Element nb_el = cc.atoms[nb.idx].el;
    if (nb_el == El::H) { has_h = true; continue; }
    if (nb_el != El::B && nb_el != El::C) return false;
    ++non_h;
  }
  return !has_h && non_h == 5;
}

} // namespace

bool has_carborane_seed(const ChemComp& cc, const AceBondAdjacency& adj) {
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    if (cc.atoms[i].el == El::H)
      continue;
    int b_count = 0;
    for (const auto& nb : adj[i])
      if (cc.atoms[nb.idx].el == El::B)
        ++b_count;
    if (b_count >= 4)
      return true;
  }
  return false;
}

bool is_carborane_mode_component(const ChemComp& cc, const AceBondAdjacency& adj) {
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    if (cc.atoms[i].el == El::H) continue;
    int b_count = 0;
    for (const auto& nb : adj[i])
      if (cc.atoms[nb.idx].el == El::B)
        ++b_count;
    if (b_count >= 4) return true;
  }
  return false;
}

void apply_carborane_mode(ChemComp& cc, bool no_angles) {
  AceGraphView initial_graph = make_ace_graph_view(cc);
  const AceBondAdjacency& initial_adj = initial_graph.adjacency;
  std::vector<size_t> h_centers;
  std::set<size_t> seen_centers;
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom), idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 >= 0 && seen_centers.insert((size_t) idx1).second && is_carborane_h_center(cc, initial_adj, (size_t) idx1))
      h_centers.push_back((size_t) idx1);
    if (idx2 >= 0 && seen_centers.insert((size_t) idx2).second && is_carborane_h_center(cc, initial_adj, (size_t) idx2))
      h_centers.push_back((size_t) idx2);
  }
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    if (seen_centers.insert(i).second && is_carborane_h_center(cc, initial_adj, i))
      h_centers.push_back(i);

  size_t h_serial = cc.atoms.size();
  for (size_t center_idx : h_centers) {
    std::string h_id = cat("H", h_serial++);
    while (cc.find_atom(h_id) != cc.atoms.end()) h_id = cat("H", h_serial++);
    cc.atoms.push_back(ChemComp::Atom{h_id, "", El::H, 0.0f, "H", "", Position()});
    cc.rt.bonds.push_back({{1, cc.atoms[center_idx].id}, {1, h_id}, BondType::Single, false, 1.10, 0.01, 1.10, 0.01});
  }

  for (auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom), idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0) continue;
    Element el1 = cc.atoms[(size_t) idx1].el, el2 = cc.atoms[(size_t) idx2].el;
    if (el1 == El::H || el2 == El::H || ((el1 == El::B || el1 == El::C) && (el2 == El::B || el2 == El::C))) {
      bond.type = BondType::Single; bond.aromatic = false;
      bond.value = bond.value_nucleus = (el1 == El::H || el2 == El::H) ? 1.10 : 1.55;
      bond.esd = bond.esd_nucleus = 0.01;
    }
  }
  cc.rt.torsions.clear(); cc.rt.chirs.clear(); cc.rt.planes.clear();
  if (no_angles) return;
  cc.rt.angles.clear();
  AceGraphView graph = make_ace_graph_view(cc);
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom), idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0) continue;
    size_t h_idx = SIZE_MAX, c_idx = SIZE_MAX;
    if (cc.atoms[(size_t) idx1].el == El::H) { h_idx = (size_t) idx1; c_idx = (size_t) idx2; }
    else if (cc.atoms[(size_t) idx2].el == El::H) { h_idx = (size_t) idx2; c_idx = (size_t) idx1; }
    else continue;
    for (const auto& nb : graph.adjacency[c_idx])
      if (nb.idx != h_idx && cc.atoms[nb.idx].el != El::H)
        cc.rt.angles.push_back({{1, cc.atoms[h_idx].id}, {1, cc.atoms[c_idx].id}, {1, cc.atoms[nb.idx].id}, 118.0, 3.0});
  }
}

void apply_mixed_carborane_mode(ChemComp& cc, bool no_angles, const std::string& tables_dir) {
  AceGraphView initial_graph = make_ace_graph_view(cc);
  std::set<size_t> cb_atoms = collect_carborane_cluster_atoms(cc, initial_graph.adjacency);
  std::set<size_t> cb_match_atoms = collect_carborone_match_atoms(cc, initial_graph.adjacency);
  if (cb_atoms.empty()) return;
  const std::set<size_t>& cb_restraint_atoms = cb_match_atoms.empty() ? cb_atoms : cb_match_atoms;

  size_t base_atom_count = cc.atoms.size();
  std::set<size_t> broken_atoms;
  std::set<std::string> reserved_tmp_h_ids;
  size_t tmp_h_serial = base_atom_count;
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom), idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0 || cc.atoms[(size_t) idx1].el == El::H || cc.atoms[(size_t) idx2].el == El::H) continue;
    if (cb_atoms.count((size_t) idx1) != cb_atoms.count((size_t) idx2)) {
      broken_atoms.insert((size_t) idx1); broken_atoms.insert((size_t) idx2);
      reserved_tmp_h_ids.insert(cat("H", ++tmp_h_serial));
    }
  }

  std::vector<size_t> h_centers;
  std::set<size_t> seen_centers;
  auto is_candidate_center = [&](size_t i) {
    if (!cb_atoms.count(i) || broken_atoms.count(i)) return false;
    if (cc.atoms[i].el != El::B && cc.atoms[i].el != El::C) return false;
    for (const auto& nb : initial_graph.adjacency[i]) if (cc.atoms[nb.idx].el == El::H) return false;
    return (int)initial_graph.adjacency[i].size() == 5;
  };
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom), idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 >= 0 && seen_centers.insert((size_t) idx1).second && is_candidate_center((size_t) idx1)) h_centers.push_back((size_t) idx1);
    if (idx2 >= 0 && seen_centers.insert((size_t) idx2).second && is_candidate_center((size_t) idx2)) h_centers.push_back((size_t) idx2);
  }
  for (size_t i = 0; i < base_atom_count; ++i) if (seen_centers.insert(i).second && is_candidate_center(i)) h_centers.push_back(i);
  if (h_centers.empty()) return;

  size_t h_serial = base_atom_count;
  std::vector<std::pair<std::string, std::string>> added_h;
  for (size_t center_idx : h_centers) {
    std::string h_id = cat("H", h_serial++);
    if (reserved_tmp_h_ids.count(h_id) || cc.find_atom(h_id) != cc.atoms.end()) h_id += "a";
    cc.atoms.push_back(ChemComp::Atom{h_id, "", El::H, 0.0f, "H", "", Position()});
    cc.rt.bonds.push_back({{1, cc.atoms[center_idx].id}, {1, h_id}, BondType::Single, false, 1.10, 0.01, 1.10, 0.01});
    added_h.emplace_back(cc.atoms[center_idx].id, h_id);
  }

  (void) apply_carborone_template_bonds(cc, cb_match_atoms, tables_dir);
  for (auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom), idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0) continue;
    bool in1 = cb_restraint_atoms.count((size_t) idx1) != 0, in2 = cb_restraint_atoms.count((size_t) idx2) != 0;
    if (!in1 && !in2) continue;
    Element el1 = cc.atoms[(size_t) idx1].el, el2 = cc.atoms[(size_t) idx2].el;
    if (el1 == El::H || el2 == El::H) {
      bond.type = BondType::Single; bond.aromatic = false;
      bond.value = bond.value_nucleus = 1.10; bond.esd = bond.esd_nucleus = 0.01;
    } else if (in1 && in2) {
      bond.type = BondType::Single; bond.aromatic = false;
      bond.value = bond.value_nucleus = (el1.is_metal() || el2.is_metal()) ? 2.10 : 1.55;
      bond.esd = bond.esd_nucleus = 0.01;
    }
  }
  std::set<std::string> h_center_ids;
  for (const auto& it : added_h) h_center_ids.insert(it.first);
  vector_remove_if(cc.rt.angles, [&](const Restraints::Angle& a) { return h_center_ids.count(a.id2.atom) != 0; });
  if (!no_angles) {
    AceGraphView graph = make_ace_graph_view(cc);
    for (const auto& it : added_h) {
      int c_idx = cc.find_atom_index(it.first), h_idx = cc.find_atom_index(it.second);
      if (c_idx < 0 || h_idx < 0) continue;
      for (const auto& nb : graph.adjacency[c_idx])
        if (nb.idx != (size_t) h_idx && cc.atoms[nb.idx].el != El::H)
          cc.rt.angles.push_back({{1, it.second}, {1, it.first}, {1, cc.atoms[nb.idx].id}, 118.0, 3.0});
    }
  }
}

} // namespace gemmi
