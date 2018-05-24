// Copyright 2017 Global Phasing Ltd.
//
// Read PDB format into a Structure from model.hpp.
// Based on the format spec:
// https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
// + support for two-character chain IDs (columns 21 and 22)
// + read segment ID (columns 73-76)
// + ignore atom serial number (compatible with the cctbx hybrid-36 extension)
// + hybrid-36 sequence id for sequences longer than 9999 (no such examples)

#ifndef GEMMI_PDB_HPP_
#define GEMMI_PDB_HPP_

#include <algorithm>  // for find_if_not, swap
#include <cassert>    // for assert
#include <cctype>     // for isspace
#include <cstdio>     // for FILE, size_t
#include <cstdlib>    // for strtol
#include <cstring>    // for memcpy, strstr
#include <map>        // for map
#include <string>     // for string
#include <vector>     // for vector

#include "model.hpp"
#include "util.hpp"
#include "fileutil.hpp" // for path_basename, file_open
#include "resinfo.hpp"
#include "polyheur.hpp" // for check_polymer_type

namespace gemmi {

namespace pdb_impl {

inline std::string rtrimmed(std::string s) {
  auto p = std::find_if_not(s.rbegin(), s.rend(),
                            [](int c) { return std::isspace(c); });
  s.erase(p.base(), s.end());
  return s;
}

inline int read_int(const char* p, int field_length) {
  int sign = 1;
  int n = 0;
  int i = 0;
  while (i < field_length && std::isspace(p[i]))
    ++i;
  if (p[i] == '-') {
    ++i;
    sign = -1;
  } else if (p[i] == '+') {
    ++i;
  }
  for (; i < field_length && p[i] >= '0' && p[i] <= '9'; ++i) {
    n = n * 10 + (p[i] - '0');
  }
  return sign * n;
}

template<int N> int read_base36(const char* p) {
  char zstr[N+1] = {0};
  std::memcpy(zstr, p, N);
  return std::strtol(zstr, nullptr, 36);
}

inline double read_double(const char* p, int field_length) {
  int sign = 1;
  double d = 0;
  int i = 0;
  while (i < field_length && std::isspace(p[i]))
    ++i;
  if (p[i] == '-') {
    ++i;
    sign = -1;
  } else if (p[i] == '+') {
    ++i;
  }
  for (; i < field_length && p[i] >= '0' && p[i] <= '9'; ++i)
    d = d * 10 + (p[i] - '0');
  if (i < field_length && p[i] == '.') {
    double mult = 0.1;
    for (++i; i < field_length && p[i] >= '0' && p[i] <= '9'; ++i, mult *= 0.1)
      d += mult * (p[i] - '0');
  }
  return sign * d;
}

inline std::string read_string(const char* p, int field_length) {
  // left trim
  while (field_length != 0 && std::isspace(*p)) {
    ++p;
    --field_length;
  }
  // EOL/EOF ends the string
  for (int i = 0; i < field_length; ++i)
    if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
      field_length = i;
      break;
    }
  // right trim
  while (field_length != 0 && std::isspace(p[field_length-1]))
    --field_length;
  return std::string(p, field_length);
}

// Compare the first 4 letters of s, ignoring case, with uppercase record.
// Both args must have at least 3+1 chars. ' ' and NUL are equivalent in s.
inline bool is_record_type(const char* s, const char* record) {
  return ((s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020) ==
          (record[0] << 24 | record[1] << 16 | record[2] << 8 | record[3]);
}

inline void setup_entities(Structure& st) {
  const static std::string WAT = "water";
  bool has_water = false;
  for (Model& model : st.models)
    for (Chain& chain : model.chains) {
      if (chain.entity_id == WAT) {
        has_water = true;
        continue;
      }
      assert(!chain.name.empty());
      if (chain.residues.size() == 1 && st.get_entity_of(chain) == nullptr)
        chain.entity_id = chain.residues[0].name + "!";
      Entity& ent = st.find_or_add_entity(chain.entity_id);
      if (chain.residues.size() == 1) {
        if (ent.entity_type == EntityType::Unknown)
          ent.entity_type = EntityType::NonPolymer;
      } else {
        PolymerType pt = check_polymer_type(chain.residues);
        if (pt == PolymerType::Unknown) {
          ent.entity_type = EntityType::NonPolymer;
        } else {
          ent.entity_type = EntityType::Polymer;
          ent.polymer_type = pt;
        }
      }
    }
  if (has_water) {
    Entity& ent = st.find_or_add_entity(WAT);
    ent.entity_type = EntityType::Water;
  }
}

inline void deduplicate_entities(Structure& st) {
  for (auto i = st.entities.begin(); i != st.entities.end(); ++i) {
    if (i->second.sequence.empty())
      continue;
    auto j = i;
    for (++j; j != st.entities.end(); )
      if (j->second.sequence == i->second.sequence) {
        for (Model& mod : st.models)
          for (Chain& ch : mod.chains)
            if (ch.entity_id == j->first)
              ch.entity_id = i->first;
        j = st.entities.erase(j);
      } else {
        ++j;
      }
  }
}

// The standard charge format is 2+, but some files have +2.
inline signed char read_charge(char digit, char sign) {
  if (sign == ' ' && digit == ' ')  // by far the most common case
    return 0;
  if (sign >= '0' && sign <= '9')
    std::swap(digit, sign);
  if (digit >= '0' && digit <= '9') {
    if (sign != '+' && sign != '-' && sign != '\0' && !std::isspace(sign))
      fail("Wrong format for charge: " +
           std::string(1, digit) + std::string(1, sign));
    return (digit - '0') * (sign == '-' ? -1 : 1);
  }
  // if we are here the field should be blank, but maybe better not to check
  return 0;
}

inline int read_matrix(Transform& t, char* line, int len) {
  if (len < 46)
    return 0;
  char n = line[5] - '0';
  if (n >= 1 && n <= 3) {
    t.mat[n-1][0] = read_double(line+10, 10);
    t.mat[n-1][1] = read_double(line+20, 10);
    t.mat[n-1][2] = read_double(line+30, 10);
    t.vec.at(n-1) = read_double(line+45, 10);
  }
  return n;
}

inline ResidueId read_res_id(const char* seq_id, const char* name) {
  ResidueId rid;
  // We support hybrid-36 extension, although it is never used in practice
  // as 9999 residues per chain are enough.
  rid.seq_num = seq_id[0] < 'A' ? read_int(seq_id, 4)
                                : read_base36<4>(seq_id) - 466560 + 10000;
  rid.icode = seq_id[4];
  rid.name = read_string(name, 3);
  return rid;
}

struct FileInput {
  std::FILE* f;
  char* gets(char* line, int size) { return std::fgets(line, size, f); }
  int getc() { return std::fgetc(f); }
};

struct MemoryInput {
  const char* start;
  const char* end;
  char* gets(char* line, int size) {
    if (start >= end)
      return nullptr;
    const char* nl = (const char*) std::memchr(start, '\n', end - start);
    size_t len = nl && nl + 1 - start < size ? nl + 1 - start : size;
    std::memcpy(line, start, len);
    start += len;
    return line;
  }
  int getc() { return start < end ? *++start : EOF; }
};

template<typename Input>
inline size_t copy_line_from_stream(char* line, int size, Input&& in) {
  if (!in.gets(line, size))
    return 0;
  size_t len = std::strlen(line);
  // If a line is longer than size we discard the rest of it.
  if (len > 0 && line[len-1] != '\n')
    for (int c = in.getc(); c != 0 && c != EOF && c != '\n'; c = in.getc())
      continue;
  return len;
}

inline SameAsu compare_link_symops(const std::string& record) {
  if (record.size() < 72)
    return SameAsu::Any;  // it could be interpreted as Yes
  if (read_string(&record[59], 6) == read_string(&record[66], 6))
    return SameAsu::Yes;
  return SameAsu::No;
}

inline
void process_conn(Structure& st, const std::vector<std::string>& conn_records) {
  int disulf_count = 0;
  int covale_count = 0;
  int metalc_count = 0;
  for (const std::string& record : conn_records) {
    if (record[0] == 'S' || record[0] == 's') { // SSBOND
      if (record.length() < 35)
        continue;
      Connection c;
      c.name = "disulf" + std::to_string(++disulf_count);
      c.type = Connection::Disulf;
      const char* r = record.c_str();
      c.atom[0].chain_name = read_string(r + 14, 2);
      c.atom[0].res_id = read_res_id(r + 17, r + 11);
      c.atom[0].atom_name = "SG";
      c.atom[1].chain_name = read_string(r + 28, 2);
      c.atom[1].res_id = read_res_id(r + 31, r + 25);
      c.atom[1].atom_name = "SG";
      c.asu = compare_link_symops(record);
      for (Model& mdl : st.models) {
        for (AtomAddress& ad : c.atom) {
          CRA cra = mdl.find_cra(ad);
          if (cra.residue) {
            // Atom name and altloc are not provided in the SSBOND record.
            // Usually it is SG (cysteine), but other disulfide bonds
            // are also possible, so if it's not CYS and SG is absent
            // we pick the first sulfur atom in the residue.
            if (!cra.atom && ad.res_id.name != "CYS")
              if (const Atom* a = cra.residue->find_by_element(El::S))
                ad.atom_name = a->name;
          }
        }
        mdl.connections.emplace_back(c);
      }
    } else if (record[0] == 'L' || record[0] == 'l') { // LINK
      if (record.length() < 57)
        continue;
      Connection c;
      // emulating names used in wwPDB mmCIFs (covaleN and metalcN)
      if (is_metal(find_element(&record[12])) ||
          is_metal(find_element(&record[42]))) {
        c.name = "metalc" + std::to_string(++metalc_count);
        c.type = Connection::MetalC;
      } else {
        c.name = "covale" + std::to_string(++covale_count);
        c.type = Connection::Covale;
      }
      for (int i : {0, 1}) {
        const char* t = record.c_str() + 30 * i;
        c.atom[i].chain_name = read_string(t + 20, 2);
        c.atom[i].res_id = read_res_id(t + 22, t + 17);
        c.atom[i].atom_name = read_string(t + 12, 4);
        c.atom[i].altloc = (t[16] == ' ' ? '\0' : t[16]);
      }
      c.asu = compare_link_symops(record);
      for (Model& mdl : st.models)
        mdl.connections.emplace_back(c);
    } else if (record[0] == 'C' || record[0] == 'c') { // CISPEP
      if (record.length() < 22)
        continue;
      const char* r = record.c_str();
      std::string cname = read_string(r + 14, 2);
      ResidueId rid = read_res_id(r + 17, r + 11);
      for (Model& model : st.models)
        for (Chain& chain : model.chains)
          if (chain.auth_name == cname)
            if (Residue* res = chain.find_residue(rid))
              res->is_cis = true;
    }
  }
}

template<typename Input>
Structure read_pdb_from_line_input(Input&& infile, const std::string& source) {
  using namespace pdb_impl;
  int line_num = 0;
  auto wrong = [&line_num](const std::string& msg) {
    fail("Problem in line " + std::to_string(line_num) + ": " + msg);
  };
  Structure st;
  st.input_format = CoorFormat::Pdb;
  st.name = path_basename(source);
  std::vector<std::string> conn_records;
  Model *model = &st.find_or_add_model("1");
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  char line[88] = {0};
  bool prev_water = false;
  Transform matrix;
  while (size_t len = copy_line_from_stream(line, 82, infile)) {
    ++line_num;
    if (is_record_type(line, "ATOM") || is_record_type(line, "HETATM")) {
      if (len < 66)
        wrong("The line is too short to be correct:\n" + std::string(line));
      std::string chain_name = read_string(line+20, 2);
      ResidueId rid = read_res_id(line+22, line+17);
      bool res_water = find_tabulated_residue(rid.name).is_water();
      if (!chain || chain_name != chain->auth_name || res_water != prev_water) {
        if (!model)
          wrong("ATOM/HETATM between models");
        if (res_water) {
          chain = &model->find_or_add_chain(chain_name + "_w");
          if (chain->entity_id.empty())
            chain->entity_id = "water";
        } else {
          int n = std::count_if(model->chains.begin(), model->chains.end(),
                  [&](const Chain& ch) { return ch.auth_name == chain_name; });
          chain = &model->add_chain(chain_name + std::to_string(n));
          chain->entity_id = (n == 0 ? chain_name : chain->name);
        }
        prev_water = res_water;
        chain->auth_name = chain_name;
        resi = nullptr;
      }
      // Non-standard but widely used 4-character segment identifier.
      // Left-justified, and may include a space in the middle.
      // The segment may be a portion of a chain or a complete chain.
      if (len > 72)
        rid.segment = read_string(line+72, 4);
      if (!resi || !resi->matches(rid)) {
        resi = chain->find_residue(rid);
        if (!resi) {
          chain->residues.emplace_back(rid);
          resi = &chain->residues.back();
        }
        resi->het_flag = line[0] & ~0x20;
      }

      Atom atom;
      atom.name = read_string(line+12, 4);
      atom.altloc = line[16] == ' ' ? '\0' : line[16];
      atom.pos.x = read_double(line+30, 8);
      atom.pos.y = read_double(line+38, 8);
      atom.pos.z = read_double(line+46, 8);
      atom.occ = (float) read_double(line+54, 6);
      atom.b_iso = (float) read_double(line+60, 6);
      bool has_elem = len > 76 && (std::isalpha(line[76]) ||
                                   std::isalpha(line[77]));
      atom.element = Element(line + (has_elem ? 76 : 12));
      atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
      resi->atoms.emplace_back(atom);

    } else if (is_record_type(line, "ANISOU")) {
      if (!model || !chain || !resi || resi->atoms.empty())
        wrong("ANISOU record not directly after ATOM/HETATM.");
      // We assume that ANISOU refers to the last atom.
      // Can it not be the case?
      Atom &atom = resi->atoms.back();
      if (atom.u11 != 0.)
        wrong("Duplicated ANISOU record or not directly after ATOM/HETATM.");
      atom.u11 = read_int(line+28, 7) * 1e-4f;
      atom.u22 = read_int(line+35, 7) * 1e-4f;
      atom.u33 = read_int(line+42, 7) * 1e-4f;
      atom.u12 = read_int(line+49, 7) * 1e-4f;
      atom.u13 = read_int(line+56, 7) * 1e-4f;
      atom.u23 = read_int(line+63, 7) * 1e-4f;

    } else if (is_record_type(line, "REMARK")) {
      st.raw_remarks.push_back(line);
      if (len > 11) {
        switch (read_int(line + 7, 3)) {
          case 2:
            if (strstr(line, "ANGSTROM"))
              st.resolution = read_double(line + 23, 7);
            break;
          default:
            // ignore all other REMARKs for now
            break;
        }
      }

    } else if (is_record_type(line, "CONECT")) {
      // ignore for now

    } else if (is_record_type(line, "SEQRES")) {
      std::string chain_name = read_string(line+10, 2);
      Entity& ent = st.find_or_add_entity(chain_name);
      ent.entity_type = EntityType::Polymer;
      for (int i = 19; i < 68; i += 4) {
        std::string res_name = read_string(line+i, 3);
        if (!res_name.empty())
          ent.sequence.emplace_back(res_name);
      }

    } else if (is_record_type(line, "HEADER")) {
      if (len > 50)
        st.info["_struct_keywords.pdbx_keywords"] =
                                    rtrimmed(std::string(line+10, 40));
      if (len > 59) { // date in PDB has format 28-MAR-07
        std::string date(line+50, 9);
        const char months[] = "JAN01FEB02MAR03APR04MAY05JUN06"
                              "JUL07AUG08SEP09OCT10NOV11DEC122222";
        const char* m = strstr(months, date.substr(3, 3).c_str());
        st.info["_pdbx_database_status.recvd_initial_deposition_date"] =
          (date[7] > '6' ? "19" : "20") + date.substr(7, 2) + "-" +
          (m ? std::string(m+3, 2) : "??") + "-" + date.substr(0, 2);
      }
      if (len > 66)
        st.info["_entry.id"] = std::string(line+62, 4);

    } else if (is_record_type(line, "TITLE")) {
      if (len > 10)
        st.info["_struct.title"] += rtrimmed(std::string(line+10, len-10-1));

    } else if (is_record_type(line, "KEYWDS")) {
      if (len > 10)
        st.info["_struct_keywords.text"] +=
                                    rtrimmed(std::string(line+10, len-10-1));

    } else if (is_record_type(line, "EXPDTA")) {
      if (len > 10)
        st.info["_exptl.method"] += rtrimmed(std::string(line+10, len-10-1));

    } else if (is_record_type(line, "CRYST1")) {
      if (len > 54)
        st.cell.set(read_double(line+6, 9),
                    read_double(line+15, 9),
                    read_double(line+24, 9),
                    read_double(line+33, 7),
                    read_double(line+40, 7),
                    read_double(line+47, 7));
      if (len > 56)
        st.spacegroup_hm = read_string(line+55, 11);
      if (len > 67) {
        std::string z = read_string(line+66, 4);
        if (!z.empty())
          st.info["_cell.Z_PDB"] = z;
      }
    } else if (is_record_type(line, "MTRIXn")) {
      if (read_matrix(matrix, line, len) == 3 && !matrix.is_identity()) {
        std::string id = read_string(line+7, 3);
        bool given = len > 59 && line[59] == '1';
        st.ncs.push_back({id, given, matrix});
        matrix.set_identity();
      }
    } else if (is_record_type(line, "MODEL")) {
      if (model && chain)
        wrong("MODEL without ENDMDL?");
      std::string name = std::to_string(read_int(line+10, 4));
      model = &st.find_or_add_model(name);
      if (!model->chains.empty())
        wrong("duplicate MODEL number: " + name);
      chain = nullptr;

    } else if (is_record_type(line, "ENDMDL")) {
      model = nullptr;
      chain = nullptr;

    } else if (is_record_type(line, "TER")) { // finishes polymer chains
      if (chain) {
        Entity& ent = st.find_or_add_entity(chain->entity_id);
        ent.entity_type = EntityType::Polymer;
      }
      chain = nullptr;

    } else if (is_record_type(line, "SCALEn")) {
      if (read_matrix(matrix, line, len) == 3) {
        st.cell.set_matrices_from_fract(matrix);
        matrix.set_identity();
      }

    } else if (is_record_type(line, "ORIGX")) {
      st.has_origx = true;
      read_matrix(st.origx, line, len);

    } else if (is_record_type(line, "SSBOND") ||
               is_record_type(line, "LINK") ||
               is_record_type(line, "CISPEP")) {
      conn_records.emplace_back(line);

    } else if (is_record_type(line, "END")) {  // NUL == ' ' & ~0x20
      break;
    }
  }

  setup_entities(st);
  deduplicate_entities(st);
  st.setup_cell_images();

  process_conn(st, conn_records);

  return st;
}

}  // namespace pdb_impl

inline Structure read_pdb_file(const std::string& path) {
  auto f = file_open(path.c_str(), "r");
  return read_pdb_from_line_input(pdb_impl::FileInput{f.get()}, path);
}

inline Structure read_pdb_from_memory(const char* data, size_t size,
                                      const std::string& name) {
  return read_pdb_from_line_input(pdb_impl::MemoryInput{data, data+size}, name);
}

inline Structure read_pdb_string(const std::string& str,
                                 const std::string& name) {
  return read_pdb_from_memory(str.c_str(), str.length(), name);
}

// A function for transparent reading of stdin and/or gzipped files.
template<typename T>
inline Structure read_pdb(T&& input) {
  if (input.is_stdin())
    return read_pdb_from_line_input(pdb_impl::FileInput{stdin}, "stdin");
  if (auto line_input = input.get_line_stream())
    return pdb_impl::read_pdb_from_line_input(line_input, input.path());
  return read_pdb_file(input.path());
}

// mmCIF files have each non-polymer residue in a separate "chain".
// This function does the same.
inline void split_nonpolymers(Structure& st) {
  for (Model& model : st.models) {
    std::vector<Chain> old_chains;
    old_chains.swap(model.chains);
    for (Chain& ch : old_chains) {
      int num = 0;
      const Entity* ent = st.get_entity(ch.entity_id);
      if (!ent || ent->entity_type != EntityType::NonPolymer ||
          ch.residues.size() < 2) {
        ch.name = ch.auth_name + std::to_string(num++);
        while (model.find_chain(ch.name) != nullptr)
          ch.name = ch.auth_name + std::to_string(num++);
        model.chains.emplace_back(std::move(ch));
      } else {
        for (Residue& res : ch.residues) {
          std::string name = ch.auth_name + std::to_string(num++);
          while (model.find_chain(name) != nullptr)
            name = ch.auth_name + std::to_string(num++);
          model.chains.emplace_back(name);
          model.chains.back().auth_name = ch.auth_name;
          model.chains.back().entity_id = res.name + "!";
          model.chains.back().residues.emplace_back(std::move(res));
        }
      }
    }
  }
  st.entities.clear();
  pdb_impl::setup_entities(st);
  pdb_impl::deduplicate_entities(st);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
