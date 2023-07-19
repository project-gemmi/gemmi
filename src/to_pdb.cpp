// Copyright 2017-2023 Global Phasing Ltd.

#include <gemmi/to_pdb.hpp>

#include <cassert>
#include <cctype>         // for isdigit
#include <cstring>        // for memset, memcpy
#include <algorithm>
#include <array>
#include <sstream>       // for ostringstream
#include <unordered_map>

#include <gemmi/fail.hpp>       // for fail
#include <gemmi/sprintf.hpp>
#include <gemmi/resinfo.hpp>
#include <gemmi/util.hpp>

namespace gemmi {

#define WRITE(...) do { \
    snprintf_z(buf, 82, __VA_ARGS__); \
    buf[80] = '\n'; \
    os.write(buf, 81); \
  } while(0)

#define WRITEU(...) do { \
    snprintf_z(buf, 82, __VA_ARGS__); \
    buf[80] = '\n'; \
    for (int i_ = 0; i_ != 80; i_++) \
      if (buf[i_] >= 'a' && buf[i_] <= 'z') buf[i_] -= 0x20; \
    os.write(buf, 81); \
  } while(0)

#define WRITELN(...) do { \
    int length__ = snprintf_z(buf, 82, __VA_ARGS__); \
    if (length__ < 80) \
      std::memset(buf + length__, ' ', 80 - length__); \
    buf[80] = '\n'; \
    os.write(buf, 81); \
  } while(0)

namespace {

bool use_hetatm(const Residue& res) {
  if (res.het_flag == 'H')
    return true;
  if (res.het_flag == 'A')
    return false;
  if (res.entity_type == EntityType::Branched ||
      res.entity_type == EntityType::NonPolymer ||
      res.entity_type == EntityType::Water)
    return true;
  return !find_tabulated_residue(res.name).is_standard();
}

// works for non-negative values only
inline void base36_encode(char* buffer, int width, int value) {
  const char base36[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  buffer[width] = '\0';
  do {
    buffer[--width] = base36[value % 36];
    value /= 36;
  } while (value != 0 && width != 0);
  while (width != 0)
    buffer[--width] = ' ';
}

// based on http://cci.lbl.gov/hybrid_36/
inline std::array<char,8> encode_serial_in_hybrid36(int serial) {
  std::array<char,8> str;
  assert(serial >= 0);
  if (serial < 100000)
    to_chars_z(str.data(), str.data() + 8, serial);
  else
    base36_encode(str.data(), 5, serial + (10 * 36 * 36 * 36 * 36 - 100000));
  return str;
}

inline std::array<char,8> write_seq_id(const SeqId& seqid) {
  std::array<char,8> str;
  char* ptr = str.data();
  if (*seqid.num > -1000 && *seqid.num < 10000) {
    ptr = to_chars_z(ptr, ptr + 5, *seqid.num);
  } else if (seqid.num.has_value()) {
    base36_encode(ptr, 4, *seqid.num + (10 * 36 * 36 * 36 - 10000));
    ptr += 4;
  }
  *ptr++ = seqid.icode;
  *ptr = '\0';
  return str;
}

inline const char* find_last_break(const char *str, int max_len) {
  int last_break = 0;
  for (int i = 0; i < max_len; i++) {
    if (str[i] == '\0')
      return str + i;
    if (str[i] == ' ' || str[i] == '-')
      last_break = i + 1;
  }
  return str + (last_break != 0 ? last_break : max_len);
}

// Write record with possible continuation lines, with the format:
// 1-6 record name, 8-10 continuation, 11-lastcol string.
inline void write_multiline(std::ostream& os, const char* record_name,
                            const std::string& text, int lastcol) {
  if (text.empty())
    return;
  char buf[88]; // a few bytes extra, just in case
  const char *start = text.c_str();
  const char *end = find_last_break(start, lastcol-10);
  WRITEU("%-6s    %-70.*s", record_name, static_cast<int>(end-start), start);
  for (int n = 2; n < 1000 && *end != '\0'; ++n) {
    start = end;
    end = find_last_break(start, lastcol-11);
    int len = int(end - start);
    WRITEU("%-6s %3d %-69.*s", record_name, n, len, start);
  }
}

inline void write_cryst1(const Structure& st, std::ostream& os) {
  char buf[88];
  const UnitCell& cell = st.cell;
  WRITE("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s          ",
        cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma,
        st.spacegroup_hm.empty() ? "P 1" : st.spacegroup_hm.c_str(),
        st.get_info("_cell.Z_PDB").c_str());
}

inline void write_ncs_op(const NcsOp& op, std::ostream& os) {
  char buf[88];
  for (int i = 0; i < 3; ++i) {
    WRITE("MTRIX%d %3.3s%10.6f%10.6f%10.6f %14.5f    %-21c", i+1,
          op.id.c_str(), op.tr.mat[i][0], op.tr.mat[i][1], op.tr.mat[i][2],
          op.tr.vec.at(i), op.given ? '1' : ' ');
  }
}

inline void write_ncs(const Structure& st, std::ostream& os) {
  if (st.ncs.empty())
    return;
  auto identity = st.info.find("_struct_ncs_oper.id");
  if (identity != st.info.end() &&
      !in_vector_f([&](const NcsOp& op) { return op.id == identity->second; }, st.ncs))
    write_ncs_op(NcsOp{identity->second, true, {}}, os);
  for (const NcsOp& op : st.ncs)
    write_ncs_op(op, os);
}

inline void write_remarks(const Structure& st, std::ostream& os) {
  char buf[88];
  if (st.resolution > 0) {
    WRITE("%-80s", "REMARK   2");
    WRITE("REMARK   2 RESOLUTION. %7.2f %-49s", st.resolution, "ANGSTROMS.");
  }
  if (!st.assemblies.empty()) {
    const char* preface[] = {
      "REMARK 350",
      "REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN",
      "REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE",
      "REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS",
      "REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND",
      "REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN."
    };
    for (const char* line : preface)
      WRITE("%-80s", line);
    int counter = 0;
    for (const Assembly& assem : st.assemblies) {
      WRITE("%-80s", "REMARK 350");
      WRITE("REMARK 350 BIOMOLECULE: %-56d", ++counter);
      if (assem.author_determined)
        WRITEU("REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: %-34s",
               assem.oligomeric_details.c_str());
      if (assem.software_determined) {
        WRITEU("REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: %-27s",
               assem.oligomeric_details.c_str());
        if (!assem.software_name.empty())
          WRITEU("REMARK 350 SOFTWARE USED: %-54s",
                 assem.software_name.c_str());
        if (!std::isnan(assem.absa))
          WRITELN("REMARK 350 TOTAL BURIED SURFACE AREA: %.0f ANGSTROM**2",
                  assem.absa);
        if (!std::isnan(assem.ssa))
          WRITELN("REMARK 350 SURFACE AREA OF THE COMPLEX: %.0f ANGSTROM**2",
                  assem.ssa);
        if (!std::isnan(assem.more))
          WRITELN("REMARK 350 CHANGE IN SOLVENT FREE ENERGY: %.1f KCAL/MOL",
                  assem.more);
      }
      int oper_cnt = 0;
      for (const Assembly::Gen& gen : assem.generators) {
        std::string chains_str;
        if (!gen.chains.empty()) {
          chains_str = join_str(gen.chains, ", ");
        } else {
          // subchains -> chains
          std::vector<std::string> chains;
          for (const Chain& ch : st.first_model().chains)
            if (!ch.residues.empty() && !in_vector(ch.name, chains) &&
                in_vector(ch.residues[0].subchain, gen.subchains))
              chains.push_back(ch.name);
          chains_str = join_str(chains, ", ");
        }
        size_t end = chains_str.length();
        if (end >= 30)
          end = chains_str.rfind(' ', 29);
        WRITELN("REMARK 350 APPLY THE FOLLOWING TO CHAINS: %s",
                chains_str.substr(0, end).c_str());
        while (end < chains_str.length()) {
          size_t begin = end + 1;
          end = chains_str.length();
          if (end - begin >= 30)
            end = chains_str.rfind(' ', begin + 29);
          WRITELN("REMARK 350                    AND CHAINS: %s",
                  chains_str.substr(begin, end - begin).c_str());
        }
        for (const Assembly::Operator& oper : gen.operators) {
          ++oper_cnt;
          const Transform& tr = oper.transform;
          for (int i = 0; i < 3; ++i)
            WRITE("REMARK 350   "
                  "BIOMT%d %3d%10.6f%10.6f%10.6f %14.5f            ",
                  i+1, oper_cnt,
                  tr.mat[i][0], tr.mat[i][1], tr.mat[i][2], tr.vec.at(i));
        }
      }
    }
  }
}

inline void write_chain_atoms(const Chain& chain, std::ostream& os,
                              int& serial, PdbWriteOptions opt) {
  char buf[88];
  buf[0] = '\0';
  if (chain.name.length() > 2)
    fail("long chain name: " + chain.name);
  for (const Residue& res : chain.residues) {
    bool as_het = use_hetatm(res);
    for (const Atom& a : res.atoms) {
      //  1- 6  6s  record name
      //  7-11  5d  integer serial
      // 12     1   -
      // 13-16  4s  atom name (from 13 only if 4-char or 2-char symbol)
      // 17     1c  altloc
      // 18-20  3s  residue name
      // 21     1   -
      // 22     1s  chain
      // 23-26  4d  integer residue sequence number
      // 27     1c  insertion code
      // 28-30  3   -
      // 31-38  8f  x (8.3)
      // 39-46  8f  y
      // 47-54  8f  z
      // 55-60  6f  occupancy (6.2)
      // 61-66  6f  temperature factor (6.2)
      // 67-76  6   -
      // 73-76      segment identifier, left-justified (non-standard)
      // 77-78  2s  element symbol, right-justified
      // 79-80  2s  charge
      WRITE("%-6s%5s %-4.4s%c%3.3s"
            "%2s%5s   %8.3f%8.3f%8.3f"
            "%6.2f%6.2f      %-4.4s%2s%c%c",
            as_het ? "HETATM" : "ATOM",
            encode_serial_in_hybrid36(++serial).data(),
            a.padded_name().c_str(),
            a.altloc ? std::toupper(a.altloc) : ' ',
            res.name.c_str(),
            chain.name.c_str(),
            write_seq_id(res.seqid).data(),
            // We want to avoid negative zero and round the numbers up
            // if they originally had one digit more and that digit was 5.
            a.pos.x > -5e-4 && a.pos.x < 0 ? 0 : a.pos.x + 1e-10,
            a.pos.y > -5e-4 && a.pos.y < 0 ? 0 : a.pos.y + 1e-10,
            a.pos.z > -5e-4 && a.pos.z < 0 ? 0 : a.pos.z + 1e-10,
            // Occupancy is stored as single prec, but we know it's <= 1,
            // so no precision is lost even if it had 6 digits after dot.
            a.occ + 1e-6,
            // B is harder to get rounded right. It is stored as float,
            // and may be given with more than single precision in mmCIF
            // If it was originally %.5f (5TIS) we need to add 0.5 * 10^-5.
            std::min(a.b_iso + 0.5e-5, 999.99),
            res.segment.c_str(),
            a.element.uname(),
            // Charge is written as 1+ or 2-, etc, or just empty space.
            // Sometimes PDB files have explicit 0s (5M05); we ignore them.
            a.charge ? a.charge > 0 ? '0'+a.charge : '0'-a.charge : ' ',
            a.charge ? a.charge > 0 ? '+' : '-' : ' ');
      ((Atom&) a).serial = serial;
      if (a.aniso.nonzero()) {
        // re-using part of the buffer
        std::memcpy(buf, "ANISOU", 6);
        const double eps = 1e-6;
        snprintf_z(buf+28, 43, "%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f",
                   a.aniso.u11*1e4 + eps, a.aniso.u22*1e4 + eps,
                   a.aniso.u33*1e4 + eps, a.aniso.u12*1e4 + eps,
                   a.aniso.u13*1e4 + eps, a.aniso.u23*1e4 + eps);
        buf[28+42] = ' ';
        buf[80] = '\n';
        os.write(buf, 81);
      }
    }
    if (opt.ter_records && buf[0] != '\0' &&
        (opt.ter_ignores_type ? &res == &chain.residues.back()
                              : (res.entity_type == EntityType::Polymer &&
                                (&res == &chain.residues.back() ||
                                 (&res + 1)->entity_type != EntityType::Polymer)))) {
      if (opt.numbered_ter) {
        // re-using part of the buffer in the middle, e.g.:
        // TER    4153      LYS B 286
        snprintf_z(buf, 82, "TER   %5s",
                   encode_serial_in_hybrid36(++serial).data());
        std::memset(buf+11, ' ', 6);
        std::memset(buf+28, ' ', 52);
        buf[80] = '\n';
        os.write(buf, 81);
      } else {
        WRITE("%-80s", "TER");
      }
    }
  }
}

inline void write_atoms(const Structure& st, std::ostream& os,
                        PdbWriteOptions opt) {
  char buf[88];
  for (const Model& model : st.models) {
    int serial = 0;
    if (st.models.size() > 1) {
      // according to the spec model name in mmCIF may not be numeric
      std::string name = model.name;
      for (char c : name)
        if (!std::isdigit(c)) {
          name = std::to_string(&model - &st.models[0] + 1);
          break;
        }
      WRITE("MODEL %8s %65s", name.c_str(), "");
    }
    for (const Chain& chain : model.chains)
      write_chain_atoms(chain, os, serial, opt);
    if (st.models.size() > 1)
      WRITE("%-80s", "ENDMDL");
    // CONECT  (note: CONECT redundancy represents bond order)
    if (opt.conect_records) {
      std::unordered_map<int, std::vector<int>> conects;
      for (const Connection& con : st.connections)
        if (con.type == Connection::Conect) {
          const_CRA cra1 = st.models[0].find_cra(con.partner1, true);
          const_CRA cra2 = st.models[0].find_cra(con.partner2, true);
          // In special cases (LINKR gap) atoms are not there.
          if (!cra1.residue || !cra2.residue)
            continue;
          conects[cra1.atom->serial].push_back(cra2.atom->serial);
          if(con.order>=2)conects[cra1.atom->serial].push_back(cra2.atom->serial);
          if(con.order>=3)conects[cra1.atom->serial].push_back(cra2.atom->serial);
        }
      for(auto ele : conects){
          
          snprintf_z(buf, 82, "CONECT%5s",
                encode_serial_in_hybrid36(ele.first).data());
          int columnOffset = 11;
          for(int c : ele.second){
            snprintf_z(buf+columnOffset, 82, "%5s",
                  encode_serial_in_hybrid36(c).data());
            columnOffset+=5;
          }
          buf[columnOffset] = '\n';
          os.write(buf, columnOffset+1);
      }
    }
  }
}

inline void write_header(const Structure& st, std::ostream& os,
                         PdbWriteOptions opt) {
  const std::string& entry_id = st.get_info("_entry.id");
  const char* entry_id_4 = entry_id.size() <= 4 ? entry_id.c_str() : "";
  char buf[88];
  { // header line
    const char* months = "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC???";
    const std::string& date =
      st.get_info("_pdbx_database_status.recvd_initial_deposition_date");
    std::string pdb_date;
    if (date.size() == 10) {
      unsigned month_idx = 10 * (date[5] - '0') + date[6] - '0' - 1;
      std::string month(months + 3 * std::min(month_idx, 13u), 3);
      pdb_date = date.substr(8, 2) + "-" + month + "-" + date.substr(2, 2);
    }
    // "classification" in PDB == _struct_keywords.pdbx_keywords in mmCIF
    const std::string& keywords = st.get_info("_struct_keywords.pdbx_keywords");
    if (!pdb_date.empty() || !keywords.empty() || !entry_id.empty())
      WRITEU("HEADER    %-40.40s%-9s   %-18.18s",
             keywords.c_str(), pdb_date.c_str(), entry_id.c_str());
  }
  write_multiline(os, "TITLE", st.get_info("_struct.title"), 80);
  write_multiline(os, "KEYWDS", st.get_info("_struct_keywords.text"), 79);
  std::string expdta = st.get_info("_exptl.method");
  if (expdta.empty())
    expdta = join_str(st.meta.experiments, "; ",
                      [](const ExperimentInfo& e) { return e.method; });
  write_multiline(os, "EXPDTA", expdta, 79);
  if (st.models.size() > 1)
    WRITE("NUMMDL    %-6zu %63s", st.models.size(), "");

  if (!st.raw_remarks.empty()) {
    for (const std::string& line : st.raw_remarks) {
      os << line;
      if (line.empty() || line.back() != '\n')
        os << '\n';
    }
  } else {
    write_remarks(st, os);
  }

  // DBREF[12], SEQRES
  if (!st.models.empty() && opt.seqres_records) {
    std::vector<const Entity*> entity_list;
    entity_list.reserve(st.models[0].chains.size());
    for (const Chain& ch : st.models[0].chains) {
      const Entity* entity = st.get_entity_of(ch.get_polymer());
      // If the input pdb file has no TER records the subchains and entities
      // are not setup automatically. In such case it is possible to call
      // setup_entities() to use heuristic to split polymer and ligands
      // and assign entities. But if it was not called, we may still find
      // the original SEQRES in the entity named after the chain name.
      if (!entity && st.input_format == CoorFormat::Pdb &&
          !ch.residues.empty() && ch.residues[0].subchain.empty()) {
        entity = st.get_entity(ch.name);
        if (entity && !entity->subchains.empty())
          entity = nullptr;
      }
      entity_list.push_back(entity);
    }
    // DBREF / DBREF1 / DBREF2
    for (size_t i = 0; i != entity_list.size(); ++i)
      if (const Entity* entity = entity_list[i]) {
        const Chain& ch = st.models[0].chains[i];
        for (const Entity::DbRef& dbref : entity->dbrefs) {
          bool short_record = *dbref.db_end.num < 100000 &&
                              dbref.accession_code.size() < 9 &&
                              dbref.id_code.size() < 13;
          SeqId begin = dbref.seq_begin;
          SeqId end = dbref.seq_end;
          if (!begin.num || !end.num)
            if (ConstResidueGroup polymer = ch.get_polymer()) {
              begin = polymer.label_seq_id_to_auth(dbref.label_seq_begin);
              end = polymer.label_seq_id_to_auth(dbref.label_seq_end);
            }
          snprintf_z(buf, 82, "DBREF  %4s%2s %5s %5s %-6s  ",
                     entry_id_4, ch.name.c_str(),
                     write_seq_id(begin).data(),
                     write_seq_id(end).data(),
                     dbref.db_name.c_str());
          if (dbref.db_name == "PDB" && dbref.id_code == entry_id) {
            // PDB uses self-reference for fragments that don't have real
            // reference. No idea why. In such case the same begin/end is used.
          } else {
            begin = dbref.db_begin;
            end = dbref.db_end;
          }
          if (short_record) {
            snprintf_z(buf+33, 82-33, "%-8s %-12s %5d%c %5d%c            \n",
                       dbref.accession_code.c_str(), dbref.id_code.c_str(),
                       *begin.num, begin.icode, *end.num, end.icode);
          } else {
            buf[5] = '1';  // -> DBREF1
            snprintf_z(buf+33, 82-33, "              %-33s\n",
                       dbref.id_code.c_str());
          }
          buf[80] = '\n';
          os.write(buf, 81);
          if (!short_record)
            WRITE("DBREF2 %4s%2s     %-22s     %10d  %10d             ",
                  entry_id_4, ch.name.c_str(),
                  dbref.accession_code.c_str(), *begin.num, *end.num);
        }
      }
    // SEQRES
    for (size_t i = 0; i != entity_list.size(); ++i)
      if (const Entity* entity = entity_list[i]) {
        const Chain& ch = st.models[0].chains[i];
        int row = 0;
        int col = 0;
        for (const std::string& monomers : entity->full_sequence) {
          if (col == 0)
            snprintf_z(buf, 82, "SEQRES%4d%2s%5zu %62s\n",
                       ++row, ch.name.c_str(),
                       entity->full_sequence.size(), "");
          size_t end = std::min(monomers.find(','), std::min(monomers.length(), (size_t)3));
          std::memcpy(buf + 18 + 4*col + 4-end, monomers.c_str(), end);
          if (++col == 13) {
            os.write(buf, 81);
            col = 0;
          }
        }
        if (col != 0) {
          buf[80] = '\n';
          os.write(buf, 81);
        }
      }
  }

  for (const ModRes& modres : st.mod_residues) {
    WRITE("MODRES %4s %3.3s%2s %5s %3s  %-41.41s  %-8.8s\n",
          entry_id_4,
          modres.res_id.name.c_str(),
          modres.chain_name.c_str(),
          write_seq_id(modres.res_id.seqid).data(),
          modres.parent_comp_id.c_str(),
          modres.details.c_str(),
          modres.mod_id.c_str());
  }

  for (const OldToNew& mapping : st.shortened_ccd_codes)
    WRITE("HETNAM     %3s %55s %-9s\n", mapping.new_.c_str(), "", mapping.old.c_str());

  if (!st.helices.empty()) {
    int counter = 0;
    for (const Helix& helix : st.helices) {
      if (++counter == 10000)
        counter = 0;
      // According to the PDB spec serial number can be from 1 to 999.
      // Here, we allow for up to 9999 helices by using columns 7 and 11.
      snprintf_z(buf, 82, "HELIX %4d%4d %3.3s%2s %5s %3.3s%2s %5s%2d %35d    \n",
            counter, counter,
            helix.start.res_id.name.c_str(), helix.start.chain_name.c_str(),
            write_seq_id(helix.start.res_id.seqid).data(),
            helix.end.res_id.name.c_str(), helix.end.chain_name.c_str(),
            write_seq_id(helix.end.res_id.seqid).data(),
            (int) helix.pdb_helix_class, helix.length);
      if (helix.length < 0) // make 72-76 blank if the length is not given
        std::memset(buf+71, ' ', 5);
      buf[80] = '\n';
      os.write(buf, 81);
    }
  }

  if (!st.sheets.empty()) {
    for (const Sheet& sheet : st.sheets) {
      int strand_counter = 0;
      for (const Sheet::Strand& strand : sheet.strands) {
        const AtomAddress& a2 = strand.hbond_atom2;
        const AtomAddress& a1 = strand.hbond_atom1;
        // H-bond atom names are expected to be O and N
        WRITE("SHEET%5d %3.3s%2zu %3.3s%2s%5s %3.3s%2s%5s%2d  %-3s%3.3s%2s%5s "
              " %-3s%3.3s%2s%5s          ",
              ++strand_counter, sheet.name.c_str(), sheet.strands.size(),
              strand.start.res_id.name.c_str(), strand.start.chain_name.c_str(),
              write_seq_id(strand.start.res_id.seqid).data(),
              strand.end.res_id.name.c_str(), strand.end.chain_name.c_str(),
              write_seq_id(strand.end.res_id.seqid).data(), strand.sense,
              a2.atom_name.c_str(), a2.res_id.name.c_str(),
              a2.chain_name.c_str(),
              a2.res_id.seqid.num ? write_seq_id(a2.res_id.seqid).data() : "",
              a1.atom_name.c_str(), a1.res_id.name.c_str(),
              a1.chain_name.c_str(),
              a1.res_id.seqid.num ? write_seq_id(a1.res_id.seqid).data() : "");
      }
    }
  }

  if (!st.models.empty()) {
    // SSBOND  (note: uses only the first model and primary conformation)
    if (opt.ssbond_records) {
      int counter = 0;
      for (const Connection& con : st.connections)
        if (con.type == Connection::Disulf) {
          const_CRA cra1 = st.models[0].find_cra(con.partner1, true);
          const_CRA cra2 = st.models[0].find_cra(con.partner2, true);
          if (!cra1.atom || !cra2.atom)
            continue;
          NearestImage im = st.cell.find_nearest_image(cra1.atom->pos,
                                                       cra2.atom->pos, con.asu);
          if (++counter == 10000)
            counter = 0;
          WRITE("SSBOND%4d %3.3s%2s %5s   %3.3s%2s %5s %28s %6s %5.2f  ",
             counter,
             cra1.residue->name.c_str(), cra1.chain->name.c_str(),
             write_seq_id(cra1.residue->seqid).data(),
             cra2.residue->name.c_str(), cra2.chain->name.c_str(),
             write_seq_id(cra2.residue->seqid).data(),
             "1555", im.symmetry_code(false).c_str(), im.dist());
        }
    }

    // LINK  (note: uses only the first model and primary conformation)
    if (opt.link_records) {
      for (const Connection& con : st.connections)
        if (con.type == Connection::Covale || con.type == Connection::MetalC ||
            con.type == Connection::Unknown) {
          const_CRA cra1 = st.models[0].find_cra(con.partner1, true);
          const_CRA cra2 = st.models[0].find_cra(con.partner2, true);
          // In special cases (LINKR gap) atoms are not there.
          if (!cra1.residue || !cra2.residue)
            continue;
          std::string im_pdb_symbol, im_dist_str;
          bool im_same_asu = true;
          if (cra1.atom && cra2.atom) {
            NearestImage im = st.cell.find_nearest_image(cra1.atom->pos,
                                                         cra2.atom->pos, con.asu);
            im_pdb_symbol = im.symmetry_code(false);
            im_dist_str = to_str_prec<2>(im.dist());
            im_same_asu = im.same_asu();
          }
          // Pdb spec: "sym1 and sym2 are right justified and are given as
          // blank when the identity operator (and no cell translation) is
          // to be applied to the atom." But all files from wwPDB have
          // 1555 not blank, so here we also write 1555,
          // except for LINKR (Refmac variant of LINK).
          snprintf_z(buf, 82, "LINK        %-4s%c%3.3s%2s%5s   "
                "            %-4s%c%3.3s%2s%5s  %6s %6s %5s  \n",
                cra1.atom ? cra1.atom->padded_name().c_str() : "",
                cra1.atom && cra1.atom->altloc ? std::toupper(cra1.atom->altloc) : ' ',
                cra1.residue->name.c_str(),
                con.partner1.chain_name.c_str(),
                write_seq_id(cra1.residue->seqid).data(),
                cra2.atom ? cra2.atom->padded_name().c_str() : "",
                cra2.atom && cra2.atom->altloc ? std::toupper(cra2.atom->altloc) : ' ',
                cra2.residue->name.c_str(),
                con.partner2.chain_name.c_str(),
                write_seq_id(cra2.residue->seqid).data(),
                "1555", im_pdb_symbol.c_str(), im_dist_str.c_str());
          if (opt.use_linkr && !con.link_id.empty()) {
            buf[4] = 'R';  // LINK -> LINKR
            if (im_same_asu)
              std::memset(buf+58, ' ', 14); // erase symmetry
            // overwrite distance with link_id
            snprintf_z(buf+72, 82-72, "%-8s\n", con.link_id.c_str());
          }
          buf[80] = '\n';
          os.write(buf, 81);
        }
    }

    // CISPEP
    if (opt.cispep_records) {
      int counter = 0;
      for (const CisPep& cispep : st.cispeps) {
        WRITE("CISPEP%4d %3.3s%2s %5s   %3.3s%2s %5s %9s %12.2f %20s",
              ++counter,
              cispep.partner_c.res_id.name.c_str(),
              cispep.partner_c.chain_name.c_str(),
              write_seq_id(cispep.partner_c.res_id.seqid).data(),
              cispep.partner_n.res_id.name.c_str(),
              cispep.partner_n.chain_name.c_str(),
              write_seq_id(cispep.partner_n.res_id.seqid).data(),
              st.models.size() > 1 ? cispep.model_str.c_str() : "0",
              std::isnan(cispep.reported_angle) ? 0. : cispep.reported_angle,
              "");
        if (counter == 9999)
          counter = 0;
      }
    }
  }

  if (opt.cryst1_record)
    write_cryst1(st, os);
  if (st.has_origx && !st.origx.is_identity()) {
    for (int i = 0; i < 3; ++i)
      WRITE("ORIGX%d %13.6f%10.6f%10.6f %14.5f %24s", i+1,
            st.origx.mat[i][0], st.origx.mat[i][1], st.origx.mat[i][2],
            st.origx.vec.at(i), "");
  }
  if (st.cell.explicit_matrices) {
    for (int i = 0; i < 3; ++i)
      // We add a small number to avoid negative 0.
      WRITE("SCALE%d %13.6f%10.6f%10.6f %14.5f %24s", i+1,
            st.cell.frac.mat[i][0] + 1e-15, st.cell.frac.mat[i][1] + 1e-15,
            st.cell.frac.mat[i][2] + 1e-15, st.cell.frac.vec.at(i) + 1e-15, "");
  }
  write_ncs(st, os);
}

inline void check_if_structure_can_be_written_as_pdb(const Structure& st) {
  for (const gemmi::Model& model : st.models)
    for (const gemmi::Chain& chain : model.chains)
      if (chain.name.size() > 2)
        gemmi::fail("chain name too long for the PDB format: " + chain.name);
}

} // anonymous namespace

std::string make_pdb_headers(const Structure& st) {
  check_if_structure_can_be_written_as_pdb(st);
  std::ostringstream os;
  write_header(st, os, PdbWriteOptions());
  return os.str();
}

void write_pdb(const Structure& st, std::ostream& os, PdbWriteOptions opt) {
  check_if_structure_can_be_written_as_pdb(st);
  write_header(st, os, opt);
  write_atoms(st, os, opt);
  char buf[88];
  WRITE("%-80s", "END");
}

void write_minimal_pdb(const Structure& st, std::ostream& os, PdbWriteOptions opt) {
  check_if_structure_can_be_written_as_pdb(st);
  write_cryst1(st, os);
  write_ncs(st, os);
  write_atoms(st, os, opt);
}

#undef WRITE
#undef WRITEU

} // namespace gemmi
