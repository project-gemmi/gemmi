// Copyright 2017 Global Phasing Ltd.
//
// Read PDB format into a Structure from model.hh.
// Based on the format spec:
// https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
// + support for two-character chain IDs (columns 21 and 22)
// + segment ID (columns 73-76)
// + the hybrid-36 extension from cctbx

#ifndef GEMMI_READ_PDB_HH_
#define GEMMI_READ_PDB_HH_

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include "model.hh"
#include <iostream> // temporary

namespace gemmi {
namespace mol {

class CstreamLineInput {
public:
  std::string source;
  CstreamLineInput(FILE* f, std::string src) : source(src), f_(f) {}

  size_t copy_line(char* line) {
    if (!fgets(line, 82, f_))
      return 0;
    size_t len = std::strlen(line);
    // We don't expect lines longer than 80 characters, but if one is found,
    // just discard the rest of the line.
    if (len > 0 && line[len-1] != '\n')
      for (int c = fgetc(f_); c != 0 && c != '\n'; c = fgetc(f_))
        continue;
    ++line_num_;
    return len;
  }
  size_t line_num() const { return line_num_; }
private:
  FILE* f_;
  size_t line_num_ = 0;
};

inline std::string rtrimmed(std::string s) {
  auto p = std::find_if_not(s.rbegin(), s.rend(),
                            [](int c) { return std::isspace(c); });
  s.erase(p.base(), s.end());
  return s;
}

inline int read_pdb_int(const char* p, int field_length) {
  int sign = 1;
  int n = 0;
  int i = 0;
  while (std::isspace(p[i]) && i < field_length)
    ++i;
  if (p[i] == '-') {
    ++i;
    sign = -1;
  } else if (p[i] == '+') {
    ++i;
  }
  for (; p[i] >= '0' && p[i] <= '9' && i < field_length; ++i) {
    n = n * 10 + (p[i] - '0');
  }
  return sign * n;
}

inline double read_pdb_number(const char* p, int field_length) {
  int sign = 1;
  double d = 0;
  int i = 0;
  while (std::isspace(p[i]) && i < field_length)
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

inline std::string read_pdb_string(const char* p, int field_length) {
  // left trim
  while (std::isspace(*p)) {
    ++p;
    --field_length;
    if (field_length == 0)
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

namespace internal {

inline Atom* find_atom_by_serial(int serial, Residue& res) {
  for (Atom& a : res.atoms)
    if ((std::intptr_t) a.parent == serial)
      return &a;
  return nullptr;
}

}


template<typename InputType>
Structure read_pdb_from_input(InputType&& in) {
  auto wrong = [&in](const std::string& msg) {
    throw std::runtime_error("Problem in line " + std::to_string(in.line_num())
                             + ": " + msg);
  };
  Structure st;
  st.name = gemmi::path_basename(in.source);
  Model *model = st.find_or_add_model("");
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  char line[88] = {0};
  size_t len;
  while ((len = in.copy_line(line))) {
    if (len < 78)
      wrong("The line is too short to be correct:\n" + std::string(line));
    if (is_record_type(line, "ATOM") || is_record_type(line, "HETATM")) {
      std::string chain_name = read_pdb_string(line+20, 2);
      if (!chain || chain_name != chain->auth_name) {
        chain = model->find_or_add_chain(chain_name);
        // if this chain was TER'ed we use a separate chain for the rest.
        if (chain->entity_type == EntityType::Polymer)
          chain = model->find_or_add_chain(chain_name + "_H");
        chain->auth_name = chain_name;
        resi = nullptr;
      }

      int seq_id = read_pdb_int(line+22, 4);
      char ins_code = line[26];
      std::string resi_name = read_pdb_string(line+17, 3);

      if (!resi || seq_id != resi->seq_id || seq_id == Residue::UnknownId ||
          resi_name != resi->name) {
        resi = chain->find_or_add_residue(seq_id, seq_id, ins_code, resi_name);
      }

      Atom atom;
      atom.name = read_pdb_string(line+12, 4);
      atom.altloc = line[16];
      atom.charge = 0;
      if (len > 78 && line[78] >= '0' && line[78] <= '9') {
        char sign = line[79];
        if (sign != '+' && sign != '-' && sign != 0 && !std::isspace(sign))
          wrong("Sign expected at position 80, got: " + std::string(1, sign));
        atom.charge = (line[78] - '0') * (sign == '-' ? -1 : 1);
      }
      atom.element = Element(line+76);
      atom.pos.x = read_pdb_number(line+30, 8);
      atom.pos.y = read_pdb_number(line+38, 8);
      atom.pos.z = read_pdb_number(line+46, 8);
      atom.occ = read_pdb_number(line+54, 6);
      atom.b_iso = read_pdb_number(line+60, 6);
      // temporarily use parent to store the serial number
      atom.parent = (Residue*) (std::intptr_t) read_pdb_int(line+6, 5);
      resi->atoms.emplace_back(atom);

    } else if (is_record_type(line, "ANISOU")) {
      if (!model)
        wrong("ANISOU between models");
      int serial = read_pdb_int(line+6, 5);
      Atom *a = nullptr;
      if (chain && resi) {
        // try the last atoms - works when ANISOU is directly after ATOM
        if ((std::intptr_t) resi->atoms.back().parent == serial) {
          a = &resi->atoms.back();
        } else {
          // search first in the last used residue
          a = internal::find_atom_by_serial(serial, *resi);
          // if not found, try the next residue
          if (!a && resi != &chain->residues.back())
            a = internal::find_atom_by_serial(serial, *++resi);
        }
      }
      if (!a)
        for (Chain& ch : model->chains) {
          for (Residue& r : ch.residues) {
            a = internal::find_atom_by_serial(serial, r);
            if (a) {
              resi = &r;
              chain = &ch;
              break;
            }
          }
          if (a)
            break;
        }
      if (!a)
        wrong("Atom serial number not found: " + std::to_string(serial));
      a->u11 = read_pdb_int(line+28, 7) * 1e-4;
      a->u22 = read_pdb_int(line+35, 7) * 1e-4;
      a->u33 = read_pdb_int(line+42, 7) * 1e-4;
      a->u12 = read_pdb_int(line+49, 7) * 1e-4;
      a->u13 = read_pdb_int(line+56, 7) * 1e-4;
      a->u23 = read_pdb_int(line+63, 7) * 1e-4;

    } else if (is_record_type(line, "REMARK")) {
    } else if (is_record_type(line, "CONECT")) {
    } else if (is_record_type(line, "HEADER")) {
      if (len > 50)
        st.info["_struct_keywords.pdbx_keywords"] =
                                    rtrimmed(std::string(line+10, 40));
      if (len > 59) { // date in PDB has format 28-MAR-07
        std::string date(line+50, 9);
        const char months[] = "JAN01FEB02MAR03APR04MAY05JUN06"
                              "JUL07AUG08SEP09OCT10NOV11DEC122222";
        const char* m = strstr(months, date.substr(3, 3).c_str());
        st.info["_database_PDB_rev.date_original"] =
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

    } else if (is_record_type(line, "CRYST1")) {
      if (len > 54) {
        st.cell.a = read_pdb_number(line+6, 9);
        st.cell.b = read_pdb_number(line+15, 9);
        st.cell.c = read_pdb_number(line+24, 9);
        st.cell.alpha = read_pdb_number(line+33, 7);
        st.cell.beta = read_pdb_number(line+40, 7);
        st.cell.gamma = read_pdb_number(line+47, 7);
        st.cell.calculate_matrices();
      }
      if (len > 56)
        st.sg_hm = read_pdb_string(line+55, 11);
      if (len > 67)
        st.info["_cell.Z_PDB"] = read_pdb_string(line+66, 4);

    } else if (is_record_type(line, "MTRIXn")) {

    } else if (is_record_type(line, "MODEL")) {
      if (!model->name.empty())
        wrong("MODEL without ENDMDL?");
      model->name = std::to_string(read_pdb_int(line+10, 4));
      if (st.find_or_add_model(model->name) != model)
        wrong("duplicate or misformatted MODEL number: " + model->name);
      chain = nullptr;

    } else if (is_record_type(line, "ENDMDL")) {
      if (model->name.empty())
        wrong("ENDMDL without MODEL?");
      model = st.find_or_add_model("");
      chain = nullptr;

    } else if (is_record_type(line, "TER")) {
      // TER record finishes polymer chains.
      if (chain)
        chain->entity_type = EntityType::Polymer;
      chain = nullptr;

    } else if (is_record_type(line, "SCALEn")) {

    } else if (is_record_type(line, "END")) {  // NUL == ' ' & ~0x20
      break;
    }
  }

  st.finish();
  return st;
}

inline Structure read_pdb_from_cstream(FILE* f, std::string source) {
  CstreamLineInput input(f, source);
  return read_pdb_from_input(input);
}

/*
inline Structure read_pdb_from_istream(std::istream &is) {
  IstreamLineInput input(is);
  return read_pdb_from_input(input);
}

inline Structure read_pdb_from_memory(const char* data, size_t size) {
  MemoryLineInput input(data, size);
  return read_pdb_from_input(input);
}
*/


inline Structure read_pdb(const std::string& path) {
  std::unique_ptr<FILE, decltype(&std::fclose)> f(std::fopen(path.c_str(), "r"),
                                                  &std::fclose);
  if (!f)
    throw std::runtime_error("Failed to open file: " + path);
  return read_pdb_from_cstream(f.get(), path);
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
