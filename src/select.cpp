// Copyright 2018 Global Phasing Ltd.

#include "gemmi/select.hpp"
#include <cstdlib>           // for strtol
#include <cctype>            // for isalpha
#include "gemmi/sprintf.hpp" // for to_str
#include "gemmi/atof.hpp"    // for fast_from_chars

namespace gemmi {

namespace {

[[noreturn]]
inline GEMMI_COLD void wrong_syntax(const std::string& cid, size_t pos,
                                    const char* info=nullptr) {
  std::string msg = "Invalid selection syntax";
  if (info)
    msg += info;
  if (pos != 0)
    cat_to(msg, " near \"", cid.substr(pos, 8), '"');
  cat_to(msg, ": ", cid);
  fail(msg);
}

inline int determine_omitted_cid_fields(const std::string& cid) {
  if (cid[0] == '/')
    return 0; // model
  if (std::isdigit(cid[0]) || cid[0] == '.' || cid[0] == '(' || cid[0] == '-')
    return 2; // residue
  size_t sep = cid.find_first_of("/([:;");
  if (sep == std::string::npos || cid[sep] == '/' || cid[sep] == ';')
    return 1; // chain
  if (cid[sep] == '(')
    return 2; // residue
  return 3;  // atom
}

inline Selection::List make_cid_list(const std::string& cid, size_t pos, size_t end,
                                     const char* disallowed_chars="-[]()!/*.:;") {
  Selection::List list;
  list.all = (cid[pos] == '*');
  list.inverted = (cid[pos] == '!');
  if (list.all || list.inverted)
    ++pos;
  list.list = cid.substr(pos, end - pos);
  // if a list have punctuation other than ',' something must be wrong
  size_t idx = list.list.find_first_of(disallowed_chars);
  if (idx != std::string::npos)
    wrong_syntax(cid, pos + idx, cat(" ('", list.list[idx], "' in a list)").c_str());
  return list;
}

inline void parse_cid_elements(const std::string& cid, size_t pos,
                               std::vector<char>& elements) {
  elements.clear();  // just in case
  if (cid[pos] == '*')
    return;
  bool inverted = false;
  if (cid[pos] == '!') {
    inverted = true;
    ++pos;
  }
  elements.resize((size_t)El::END, char(inverted));
  for (;;) {
    size_t sep = cid.find_first_of(",]", pos);
    if (sep == pos || sep > pos + 2) {
      if (sep == pos + 6 && cid.compare(pos, 6, "metals", 6) == 0) {
        for (size_t i = 0; i < elements.size(); ++i)
          if (is_metal(static_cast<El>(i)))
            elements[i] = char(!inverted);
      } else if (sep == pos + 9 && cid.compare(pos, 9, "nonmetals", 9) == 0) {
        for (size_t i = 0; i < elements.size(); ++i)
          if (!is_metal(static_cast<El>(i)))
            elements[i] = char(!inverted);
      } else {
        wrong_syntax(cid, 0, " in [...]");
      }
    } else {
      char elem_str[2] = {cid[pos], sep > pos+1 ? cid[pos+1] : '\0'};
      Element el = find_element(elem_str);
      if (el == El::X && (alpha_up(elem_str[0]) != 'X' || elem_str[1] != '\0'))
        wrong_syntax(cid, 0, " (invalid element in [...])");
      elements[el.ordinal()] = char(!inverted);
    }
    if (cid[sep] == ']')
      break;
    pos = sep + 1;
  }
}

inline Selection::SequenceId parse_cid_seqid(const std::string& cid, size_t& pos,
                                             int default_seqnum) {
  size_t initial_pos = pos;
  int seqnum = default_seqnum;
  char icode = ' ';
  if (cid[pos] == '*') {
    ++pos;
    icode = '*';
  } else if (std::isdigit(cid[pos]) || cid[pos] == '-') {
    char* endptr;
    seqnum = std::strtol(&cid[pos], &endptr, 10);
    pos = endptr - &cid[0];
  }
  if (cid[pos] == '.')
    ++pos;
  if (initial_pos != pos && (std::isalpha(cid[pos]) || cid[pos] == '*'))
    icode = cid[pos++];
  return {seqnum, icode};
}

inline Selection::AtomInequality parse_atom_inequality(const std::string& cid,
                                                       size_t pos, size_t end) {
  Selection::AtomInequality r;
  if (cid[pos] != 'q' && cid[pos] != 'b')
    wrong_syntax(cid, pos);
  r.property = cid[pos];
  ++pos;
  while (cid[pos] == ' ')
    ++pos;
  if (cid[pos] == '<')
    r.relation = -1;
  else if (cid[pos] == '>')
    r.relation = 1;
  else if (cid[pos] == '=')
    r.relation = 0;
  else
    wrong_syntax(cid, pos);
  ++pos;
  auto result = fast_from_chars(cid.c_str() + pos, r.value);
  if (result.ec != std::errc())
    wrong_syntax(cid, pos, " (expected number)");
  pos = size_t(result.ptr - cid.c_str());
  while (cid[pos] == ' ')
    ++pos;
  if (pos != end)
    wrong_syntax(cid, pos);
  return r;
}

inline bool has_inequality(const std::string& cid, size_t start, size_t end) {
  for (size_t i = start; i < end; ++i)
    if (cid[i] == '<' || cid[i] == '=' || cid[i] == '>')
      return true;
  return false;
}

inline void parse_cid(const std::string& cid, Selection& sel) {
  if (cid.empty() || (cid.size() == 1 && cid[0] == '*'))
    return;
  int omit = determine_omitted_cid_fields(cid);
  size_t sep = 0;
  size_t semi = cid.find(';');
  // model
  if (omit == 0) {
    sep = std::min(cid.find('/', 1), semi);
    if (sep != 1 && cid[1] != '*') {
      char* endptr;
      sel.mdl = std::strtol(&cid[1], &endptr, 10);
      size_t end_pos = endptr - &cid[0];
      if (end_pos != sep && end_pos != cid.size())
        wrong_syntax(cid, 0, " (at model number)");
    }
  }

  // chain
  if (omit <= 1 && sep < semi) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    sep = std::min(cid.find('/', pos), semi);
    // These characters are not really disallowed, but are unexpected.
    // "-" is expected, it's in chain IDs in bioassembly files from RCSB.
    const char* disallowed_chars = "[]()!/*.:;";
    sel.chain_ids = make_cid_list(cid, pos, sep, disallowed_chars);
  }

  // residue; MMDB CID syntax: s1.i1-s2.i2 or *(res).ic
  // In gemmi both 14.a and 14a are accepted.
  // *(ALA). and *(ALA) and (ALA). can be used instead of (ALA) for
  // compatibility with MMDB.
  if (omit <= 2 && sep < semi) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    if (cid[pos] != '(')
      sel.from_seqid = parse_cid_seqid(cid, pos, INT_MIN);
    if (cid[pos] == '(') {
      ++pos;
      size_t right_br = cid.find(')', pos);
      sel.residue_names = make_cid_list(cid, pos, right_br);
      pos = right_br + 1;
    }
    // allow "(RES)." and "(RES).*" and "(RES)*"
    if (cid[pos] == '.')
      ++pos;
    if (cid[pos] == '*')
      ++pos;
    if (cid[pos] == '-') {
      ++pos;
      sel.to_seqid = parse_cid_seqid(cid, pos, INT_MAX);
    } else if (sel.from_seqid.seqnum != INT_MIN) {
      sel.to_seqid = sel.from_seqid;
    }
    sep = pos;
    if (cid[sep] != '/' && cid[sep] != ';' && cid[sep] != '\0')
      wrong_syntax(cid, 0, " (at residue)");
  }

  // atom;  at[el]:aloc
  if (sep < std::min(cid.size(), semi)) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    size_t end = cid.find_first_of("[:;", pos);
    if (end != pos) {
      sel.atom_names = make_cid_list(cid, pos, end);
      // Chain name can be empty, but not atom name,
      // so we interpret empty atom name as *.
      if (!sel.atom_names.inverted && sel.atom_names.list.empty())
        sel.atom_names.all = true;
      if (end == std::string::npos)
        return;
    }
    if (cid[end] == '[') {
      pos = end + 1;
      end = cid.find(']', pos);
      if (end == std::string::npos)
        wrong_syntax(cid, 0, " (no matching ']')");
      parse_cid_elements(cid, pos, sel.elements);
      ++end;
    }
    if (cid[end] == ':') {
      pos = end + 1;
      sel.altlocs = make_cid_list(cid, pos, semi);
    }
  }

  // extensions after semicolon(s)
  while (semi < cid.size()) {
    size_t pos = semi + 1;
    while (cid[pos] == ' ')
      ++pos;
    semi = std::min(cid.find(';', pos), cid.size());
    size_t end = semi;
    while (end > pos && cid[end-1] == ' ')
      --end;
    if (has_inequality(cid, pos, end)) {
      sel.atom_inequalities.push_back(parse_atom_inequality(cid, pos, end));
    } else {
      sel.entity_types = make_cid_list(cid, pos, end);
      bool inv = sel.entity_types.inverted;
      std::fill(sel.et_flags.begin(), sel.et_flags.end(), char(inv));
      for (const std::string& item : split_str(sel.entity_types.list, ',')) {
        EntityType et = EntityType::Unknown;
        if (item == "polymer")
          et = EntityType::Polymer;
        else if (item == "solvent")
          et = EntityType::Water;
        else
          wrong_syntax(cid, 0, (" at " + item).c_str());
        sel.et_flags[(int)et] = char(!inv);
      }
    }
  }
}

} // anonymous namespace

Selection::Selection(const std::string& cid) {
  parse_cid(cid, *this);
}

std::string Selection::SequenceId::str() const {
  std::string s;
  if (!empty()) {
    s = std::to_string(seqnum);
    if (icode != '*') {
      s += '.';
      if (icode != ' ')
        s += icode;
    }
  }
  return s;
}

std::string Selection::AtomInequality::str() const {
  std::string r = ";";
  r += property;
  r += relation == 0 ? '=' : relation < 0 ? '<' : '>';
  r += to_str(value);
  return r;
}

std::string Selection::str() const {
  std::string cid = "/";
  if (mdl != 0)
    cid += std::to_string(mdl);
  cid += '/';
  cid += chain_ids.str();
  cid += '/';
  cid += from_seqid.str();
  if (!residue_names.all) {
    cid += '(';
    cid += residue_names.str();
    cid += ')';
  }
  if ((!from_seqid.empty() || !to_seqid.empty()) &&
      (from_seqid.seqnum != to_seqid.seqnum || from_seqid.icode != to_seqid.icode)) {
    cid += '-';
    cid += to_seqid.str();
  }
  cid += '/';
  if (!atom_names.all)
    cid += atom_names.str();
  if (!elements.empty()) {
    cid += '[';
    bool inv = (std::count(elements.begin(), elements.end(), 1) > 64);
    if (inv)
      cid += '!';
    for (size_t i = 0; i < elements.size(); ++i)
      if (elements[i] != char(inv)) {
        cid += element_name(static_cast<El>(i));
        cid += ',';
      }
    cid.back() = ']';
  }
  if (!altlocs.all) {
    cid += ':';
    cid += altlocs.str();
  }
  if (!entity_types.all) {
    cid += ';';
    cid += entity_types.str();
  }
  for (const AtomInequality& ai : atom_inequalities)
    cid += ai.str();
  return cid;
}

} // namespace gemmi
