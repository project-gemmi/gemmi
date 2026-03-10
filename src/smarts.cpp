// Copyright 2026 Global Phasing Ltd.
//
// Lightweight SMARTS-subset matcher for ChemComp.

#include "gemmi/smarts.hpp"
#include "gemmi/ace_graph.hpp"
#include <algorithm>
#include <cctype>
#include <functional>

namespace gemmi {

namespace {

enum class SmartsBond {
  Implicit,  // default bond: single for aliphatic, aromatic for aromatic atoms
  Any,
  Single,
  Double
};

struct SmartsNode {
  bool wildcard = false;
  int aromatic = -1;  // -1: don't care, 0: aliphatic, 1: aromatic
  El element = El::X;
  int h_count = -1;   // -1 = not specified
  int degree = -1;     // -1 = not specified (total connections incl. H)
};

struct SmartsEdge {
  int a = -1;
  int b = -1;
  SmartsBond bond = SmartsBond::Single;
};

struct SmartsPattern {
  std::vector<SmartsNode> nodes;
  std::vector<SmartsEdge> edges;
  std::vector<std::vector<std::pair<int, SmartsBond>>> adj;
};

bool parse_smarts_atom(const std::string& s, size_t& pos, SmartsNode& out) {
  if (pos >= s.size())
    return false;
  if (s[pos] == '*') {
    out.wildcard = true;
    out.aromatic = -1;
    ++pos;
    return true;
  }
  auto set_symbol = [&](const std::string& symbol, bool aromatic) {
    out.wildcard = false;
    out.aromatic = aromatic ? 1 : 0;
    out.element = Element(symbol).elem;
    return out.element != El::X;
  };
  if (s[pos] == '[') {
    size_t end = s.find(']', pos + 1);
    if (end == std::string::npos)
      return false;
    std::string body = s.substr(pos + 1, end - pos - 1);
    pos = end + 1;
    size_t i = 0;
    while (i < body.size() && !std::isalpha(static_cast<unsigned char>(body[i])))
      ++i;
    if (i >= body.size())
      return false;
    bool aromatic = std::islower(static_cast<unsigned char>(body[i]));
    std::string sym(1, static_cast<char>(std::toupper(static_cast<unsigned char>(body[i]))));
    ++i;
    if (i < body.size() && std::islower(static_cast<unsigned char>(body[i])) && !aromatic) {
      sym.push_back(body[i]);
      ++i;
    }
    if (!set_symbol(sym, aromatic))
      return false;
    while (i < body.size()) {
      char bc = body[i];
      if (bc == 'H') {
        ++i;
        if (i < body.size() && std::isdigit(static_cast<unsigned char>(body[i]))) {
          out.h_count = body[i] - '0';
          ++i;
        } else {
          out.h_count = 1;
        }
      } else if (bc == 'X') {
        ++i;
        if (i < body.size() && std::isdigit(static_cast<unsigned char>(body[i]))) {
          out.degree = body[i] - '0';
          ++i;
        }
      } else {
        ++i;
      }
    }
    return true;
  }
  if (std::isalpha(static_cast<unsigned char>(s[pos]))) {
    bool aromatic = std::islower(static_cast<unsigned char>(s[pos]));
    std::string sym(1, static_cast<char>(std::toupper(static_cast<unsigned char>(s[pos]))));
    ++pos;
    if (!aromatic && pos < s.size() &&
        std::islower(static_cast<unsigned char>(s[pos])) &&
        s[pos] != 'c' && s[pos] != 'n' && s[pos] != 'o' && s[pos] != 's' &&
        s[pos] != 'p') {
      sym.push_back(s[pos]);
      ++pos;
    }
    return set_symbol(sym, aromatic);
  }
  return false;
}

bool parse_smarts_subset(const std::string& s, SmartsPattern& out) {
  out = SmartsPattern{};
  std::vector<int> branch_stack;
  int current = -1;
  SmartsBond pending = SmartsBond::Implicit;
  size_t pos = 0;
  while (pos < s.size()) {
    char c = s[pos];
    if (std::isspace(static_cast<unsigned char>(c))) {
      ++pos;
      continue;
    }
    if (c == '(') {
      if (current < 0) return false;
      branch_stack.push_back(current);
      ++pos;
      continue;
    }
    if (c == ')') {
      if (branch_stack.empty()) return false;
      current = branch_stack.back();
      branch_stack.pop_back();
      ++pos;
      continue;
    }
    if (c == '=') {
      pending = SmartsBond::Double;
      ++pos;
      continue;
    }
    if (c == '-') {
      pending = SmartsBond::Single;
      ++pos;
      continue;
    }
    if (c == '~') {
      pending = SmartsBond::Any;
      ++pos;
      continue;
    }
    SmartsNode node;
    size_t before = pos;
    if (!parse_smarts_atom(s, pos, node))
      return false;
    if (pos == before)
      return false;
    int node_id = static_cast<int>(out.nodes.size());
    out.nodes.push_back(node);
    if (current >= 0)
      out.edges.push_back({current, node_id, pending});
    current = node_id;
    pending = SmartsBond::Implicit;
  }
  if (!branch_stack.empty())
    return false;
  out.adj.assign(out.nodes.size(), {});
  for (const SmartsEdge& e : out.edges) {
    out.adj[e.a].push_back({e.b, e.bond});
    out.adj[e.b].push_back({e.a, e.bond});
  }
  return !out.nodes.empty();
}

int count_h_neighbors(const ChemComp& cc, const AceBondAdjacency& adj, size_t idx) {
  int n = 0;
  for (const AceBondNeighbor& nb : adj[idx])
    if (cc.atoms[nb.idx].el == El::H)
      ++n;
  return n;
}

bool smarts_node_matches(const SmartsNode& p, const ChemComp& cc,
                         const AceBondAdjacency& adj, size_t atom_idx) {
  const ChemComp::Atom& a = cc.atoms[atom_idx];
  if (!p.wildcard && a.el != p.element)
    return false;
  if (p.aromatic != -1) {
    bool atom_is_arom = atom_has_aromatic_bond(adj, atom_idx);
    if ((p.aromatic == 1) != atom_is_arom)
      return false;
  }
  if (p.h_count >= 0 && count_h_neighbors(cc, adj, atom_idx) != p.h_count)
    return false;
  if (p.degree >= 0 && static_cast<int>(adj[atom_idx].size()) != p.degree)
    return false;
  return true;
}

bool smarts_bond_matches(SmartsBond p, BondType bt,
                         bool both_aromatic) {
  if (p == SmartsBond::Any)
    return true;
  if (p == SmartsBond::Single)
    return bt == BondType::Single;
  if (p == SmartsBond::Double)
    return bt == BondType::Double || bt == BondType::Deloc;
  // Implicit: aromatic bond if both atoms are aromatic, otherwise single
  if (both_aromatic)
    return is_aromatic_or_deloc(bt);
  return bt == BondType::Single;
}

BondType find_cc_bond_type(const AceBondAdjacency& adj, size_t a, size_t b) {
  for (const AceBondNeighbor& nb : adj[a])
    if (nb.idx == b)
      return nb.type;
  return BondType::Unspec;
}

} // namespace

std::vector<SmartsMatch> match_smarts(const ChemComp& cc, const std::string& pattern_text) {
  SmartsPattern p;
  if (!parse_smarts_subset(pattern_text, p))
    return {};
  AceGraphView gv = make_ace_graph_view(cc);
  const size_t n = p.nodes.size();
  std::vector<SmartsMatch> results;
  std::vector<int> map_p2c(n, -1);
  std::vector<char> used(cc.atoms.size(), 0);

  std::vector<int> order(n);
  for (size_t i = 0; i < n; ++i)
    order[i] = static_cast<int>(i);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    const SmartsNode& na = p.nodes[a];
    const SmartsNode& nb = p.nodes[b];
    int sa = (na.wildcard ? 0 : 4) + (na.aromatic != -1 ? 2 : 0) +
             static_cast<int>(p.adj[a].size());
    int sb = (nb.wildcard ? 0 : 4) + (nb.aromatic != -1 ? 2 : 0) +
             static_cast<int>(p.adj[b].size());
    return sa > sb;
  });

  std::function<void(size_t)> dfs = [&](size_t depth) {
    if (depth == n) {
      results.push_back(map_p2c);
      return;
    }
    int pi = order[depth];
    for (size_t ci = 0; ci < cc.atoms.size(); ++ci) {
      if (used[ci])
        continue;
      if (!smarts_node_matches(p.nodes[pi], cc, gv.adjacency, ci))
        continue;
      bool ok = true;
      for (const auto& nb : p.adj[pi]) {
        int pj = nb.first;
        int cj = map_p2c[pj];
        if (cj < 0)
          continue;
        BondType bt = find_cc_bond_type(gv.adjacency, ci, static_cast<size_t>(cj));
        bool both_arom = p.nodes[pi].aromatic == 1 && p.nodes[pj].aromatic == 1;
        if (bt == BondType::Unspec || !smarts_bond_matches(nb.second, bt, both_arom)) {
          ok = false;
          break;
        }
      }
      if (!ok)
        continue;
      map_p2c[pi] = static_cast<int>(ci);
      used[ci] = 1;
      dfs(depth + 1);
      used[ci] = 0;
      map_p2c[pi] = -1;
    }
  };

  dfs(0);
  return results;
}

} // namespace gemmi
