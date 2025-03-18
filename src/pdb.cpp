// Copyright 2019 Global Phasing Ltd.

#include "gemmi/pdb.hpp"
#include <cctype>             // for isalpha
#include <cstdlib>            // for atoi, strtol
#include <cstring>            // for memcpy, strstr, strchr, strcmp
#include <algorithm>          // for min, swap
#include <stdexcept>          // for invalid_argument
#include <unordered_map>
#include "gemmi/atof.hpp"     // for fast_from_chars
#include "gemmi/atox.hpp"     // for is_space, is_digit
#include "gemmi/input.hpp"
#include "gemmi/metadata.hpp" // for Metadata
#include "gemmi/model.hpp"    // for Structure, impl::find_or_add
#include "gemmi/polyheur.hpp" // for assign_subchains
#include "gemmi/util.hpp"     // for trim_str, alpha_up, istarts_with

namespace gemmi {

namespace {

int read_int(const char* p, int field_length) {
  return string_to_int(p, false, field_length);
}

double read_double(const char* p, int field_length) {
  double d = 0.;
  // we don't check for errors here
  fast_from_chars(p, p + field_length, d);
  return d;
}

std::string read_string(const char* p, int field_length) {
  // left trim
  while (field_length != 0 && is_space(*p)) {
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
  while (field_length != 0 && is_space(p[field_length-1]))
    --field_length;
  return std::string(p, field_length);
}

template<int N> int read_base36(const char* p) {
  char zstr[N+1] = {0};
  std::memcpy(zstr, p, N);
  return std::strtol(zstr, nullptr, 36);
}

// The standard charge format is 2+, but some files have +2.
signed char read_charge(char digit, char sign) {
  if (sign == ' ' && digit == ' ')  // by far the most common case
    return 0;
  if (sign >= '0' && sign <= '9')
    std::swap(digit, sign);
  if (digit >= '0' && digit <= '9') {
    if (sign != '+' && sign != '-' && sign != '\0' && !is_space(sign))
      fail("Wrong format for charge: " +
           std::string(1, digit) + std::string(1, sign));
    return (digit - '0') * (sign == '-' ? -1 : 1);
  }
  // if we are here the field should be blank, but maybe better not to check
  return 0;
}

int read_matrix(Transform& t, const char* line, size_t len) {
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

SeqId read_seq_id(const char* str) {
  SeqId seqid;
  if (str[4] != '\r' && str[4] != '\n')
    seqid.icode = str[4];
  // We support hybrid-36 extension, although it is never used in practice
  // as 9999 residues per chain are enough.
  if (str[0] < 'A') {
    for (int i = 4; i != 0; --i, ++str)
      if (!is_space(*str)) {
        seqid.num = read_int(str, i);
        break;
      }
  } else {
    seqid.num = read_base36<4>(str) - 466560 + 10000;
  }
  return seqid;
}

ResidueId read_res_id(const char* seq_id, const char* name) {
  return {read_seq_id(seq_id), {}, read_string(name, 3)};
}

char read_altloc(char c) { return c == ' ' ? '\0' : c; }

int read_serial(const char* ptr) {
  return ptr[0] < 'A' ? read_int(ptr, 5)
                      : read_base36<5>(ptr) - 16796160 + 100000;
}

El infer_element_from_padded_name(const char* name) {
  // Old versions of the PDB format had hydrogen names such as "1HB ".
  // Some MD files use similar names for other elements ("1C4A" -> C).
  if (name[0] == ' ' || is_digit(name[0]))
    return impl::find_single_letter_element(name[1]);
  // ... or it can be "C210"
  if (is_digit(name[1]))
    return impl::find_single_letter_element(name[0]);
  if (name[3] != ' ') {
    // Atom names HXXX are ambiguous, but Hg, He, Hf, Ho and Hs (almost)
    // never have 4-character names, so H is assumed.
    if (alpha_up(name[0]) == 'H')
      return El::H;
    // Similarly Deuterium (DXXX), but here alternatives are Dy, Db and Ds.
    // Only Dysprosium is present in the PDB - in a single entry as of 2022.
    if (alpha_up(name[0]) == 'D')
      return El::D;
    // Don't try harder for now. We don't recognize names such as CG11 as C
    // (which we could; there is no Cg in the periodic table), but
    // a name such as CL20 can be either Cl (in WGW) or C (in WQH) ¯\_(ツ)_/¯
  }
  return find_element(name);
}

bool element_from_padded_name_is_ambiguous(const char* name) {
  return name[0] != ' ' && name[3] != ' ' && !is_digit(name[0]) && !is_digit(name[1]);
}

// "28-MAR-07" -> "2007-03-28"
// (we also accept less standard format "28-Mar-2007" as used by BUSTER)
// We do not check if the date is correct.
// The returned value is one of:
//   DDDD-DD-DD - possibly correct date,
//   DDDD-xx-DD - unrecognized month,
//   empty string - the digits were not there.
std::string pdb_date_format_to_iso(const std::string& date) {
  const char months[] = "JAN01FEB02MAR03APR04MAY05JUN06"
                        "JUL07AUG08SEP09OCT10NOV11DEC122222";
  if (date.size() < 9 || !is_digit(date[0]) || !is_digit(date[1]) ||
                         !is_digit(date[7]) || !is_digit(date[8]))
    return std::string();
  std::string iso = "xxxx-xx-xx";
  if (date.size() >= 11 && is_digit(date[9]) && is_digit(date[10])) {
    std::memcpy(&iso[0], &date[7], 4);
  } else {
    std::memcpy(&iso[0], (date[7] > '6' ? "19" : "20"), 2);
    std::memcpy(&iso[2], &date[7], 2);
  }
  char month[4] = {alpha_up(date[3]), alpha_up(date[4]), alpha_up(date[5]), '\0'};
  if (const char* m = std::strstr(months, month))
    std::memcpy(&iso[5], m + 3, 2);
  std::memcpy(&iso[8], &date[0], 2);
  return iso;
}


bool is_double(const char* p) {
  while (is_space(*p)) ++p;
  if (*p == '-' || *p == '+') ++p;
  while (is_digit(*p)) ++p;
  if (*p == '.') {
    ++p;
    while (is_digit(*++p)) ++p;
  }
  while (is_space(*p)) ++p;
  return *p == '\0';
}

template<size_t N>
bool same_str(const std::string& s, const char (&literal)[N]) {
  return s.size() == N - 1 && std::strcmp(s.c_str(), literal) == 0;
}

bool is_tls_item(const std::string& key) {
  return key.size() == 3 &&
    (key[0] == 'T' || key[0] == 'L' || key[0] == 'S') &&
    (key[1] == '1' || key[1] == '2' || key[1] == '3') &&
    (key[2] == '1' || key[2] == '2' || key[2] == '3');
}

// Usually we have one program per line:
//   XDS
//   XDS VERSION NOVEMBER 3, 2014
//   AIMLESS 0.5.17
// but it can also be a list of programs:
//   autoPROC (Version 1.3.0), AIMLESS, STARANISO
//   autoPROC, XDS (VERSION Jan 26, 2018)
// We assume that:
// - the name has only one word (apologies to Queen of Spades,
//   Force Field X, APEX 2 and Insight II).
// - comma not followed by a digit separates programs
// - brackets and the word VERSION are to be removed from version
// Additionally, if version has format: "something (DATE)" where
// the DATE format is either 28-MAR-07 or 28-Mar-2007, then DATE
// is put into _software.date.
void add_software(Metadata& meta, SoftwareItem::Classification type, const std::string& name) {
  for (size_t start = 0, end = 0; end != std::string::npos; start = end + 1) {
    end = name.find(',', start);
    while (end != std::string::npos &&
           name[end+1] == ' ' && is_digit(name[end+2]))
      end = name.find(',', end + 1);
    meta.software.emplace_back();
    SoftwareItem& item = meta.software.back();
    item.name = trim_str(name.substr(start, end - start));
    size_t sep = item.name.find(' ');
    if (sep != std::string::npos) {
      size_t ver_start = item.name.find_first_not_of(" (", sep + 1);
      item.version = item.name.substr(ver_start);
      item.name.resize(sep);
      if (!item.version.empty() && item.version.back() == ')') {
        size_t open_br = item.version.find('(');
        if (open_br == std::string::npos) {
          item.version.pop_back();
        } else if (open_br + 11 == item.version.size() ||
                   open_br + 13 == item.version.size()) {
          item.date = pdb_date_format_to_iso(item.version.substr(open_br + 1));
          if (item.date.size() == 10 && item.date[5] != 'x') {
            size_t last = item.version.find_last_not_of(' ', open_br - 1);
            item.version.resize(last + 1);
          } else {
            item.date.clear();
          }
        }
      }
      if (istarts_with(item.version, "version "))
        item.version.erase(0, 8);
    }
    item.classification = type;
  }
}

// REMARK   3   TERM                          COUNT    WEIGHT   FUNCTION.
// REMARK   3    BOND LENGTHS              : 5760   ; 2.000  ; HARMONIC
void add_restraint_count_weight(RefinementInfo& ref_info, const char* key, const char* value) {
  if (*value == 'N') // NULL instead of number
    return;
  ref_info.restr_stats.emplace_back(key);
  RefinementInfo::Restr& restr = ref_info.restr_stats.back();
  const char* endptr;
  restr.count = no_sign_atoi(value, &endptr);
  if (const char* sep = std::strchr(endptr, ';'))
    restr.weight = fast_atof(sep + 1, &endptr);
  if (const char* sep = std::strchr(endptr, ';'))
    restr.function = read_string(sep+1, 50);
}

void read_remark3_line(const char* line, Metadata& meta,
                       std::string*& possibly_unfinished_remark3) {
  // Based on:
  // www.wwpdb.org/documentation/file-format-content/format23/remark3.html
  // and analysis of PDB files.
  // In special cases, such as joint X-ray and neutron refinement 5MOO,
  // PDB file can have two REMARK 3 blocks.
  // Generally, after "REMARK   3" we have either a header-like sentence
  // or a key:value pair with a colon, or a continuation of text from the
  // previous line.
  const char* key_start = skip_blank(line + 10);
  const char* colon = std::strchr(key_start, ':');
  const char* key_end = rtrim_cstr(key_start, colon);
  std::string key(key_start, key_end);

  // multi-line continuation requires special handling
  if (possibly_unfinished_remark3) {
    if (key_start > line + 17) {
      *possibly_unfinished_remark3 += ' ';
      possibly_unfinished_remark3->append(key);
      return;
    }
    possibly_unfinished_remark3 = nullptr;
  }

  if (colon) {
    const char* value = skip_blank(colon + 1);
    const char* end = rtrim_cstr(value);
    if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
      return;
    if (same_str(key, "PROGRAM"))
      add_software(meta, SoftwareItem::Refinement, std::string(value, end));
    if (meta.refinement.empty())
      return;
    RefinementInfo& ref_info = meta.refinement.back();
    if (same_str(key, "RESOLUTION RANGE HIGH (ANGSTROMS)")) {
      ref_info.resolution_high = fast_atof(value);
    } else if (same_str(key, "RESOLUTION RANGE LOW  (ANGSTROMS)")) {
      ref_info.resolution_low = fast_atof(value);
    } else if (same_str(key, "COMPLETENESS FOR RANGE        (%)")) {
      ref_info.completeness = fast_atof(value);
    } else if (same_str(key, "NUMBER OF REFLECTIONS")) {
      ref_info.reflection_count = std::atoi(value);
    } else if (same_str(key, "CROSS-VALIDATION METHOD")) {
      ref_info.cross_validation_method = std::string(value, end);
    } else if (same_str(key, "FREE R VALUE TEST SET SELECTION")) {
      ref_info.rfree_selection_method = std::string(value, end);
    } else if (same_str(key, "R VALUE     (WORKING + TEST SET)")) {
      ref_info.r_all = fast_atof(value);
    } else if (same_str(key, "R VALUE            (WORKING SET)")) {
      ref_info.r_work = fast_atof(value);
    } else if (same_str(key, "FREE R VALUE")) {
      ref_info.r_free = fast_atof(value);
    } else if (same_str(key, "FREE R VALUE TEST SET COUNT")) {
      ref_info.rfree_set_count = atoi(value);
    } else if (same_str(key, "TOTAL NUMBER OF BINS USED")) {
      ref_info.bin_count = std::atoi(value);
    } else if (same_str(key, "BIN RESOLUTION RANGE HIGH       (A)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().resolution_high = fast_atof(value);
    } else if (same_str(key, "BIN RESOLUTION RANGE LOW        (A)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().resolution_low = fast_atof(value);
    } else if (same_str(key, "BIN COMPLETENESS (WORKING+TEST) (%)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().completeness = fast_atof(value);
    } else if (same_str(key, "REFLECTIONS IN BIN   (WORKING+TEST)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().reflection_count = std::atoi(value);
    } else if (same_str(key, "BIN R VALUE          (WORKING+TEST)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_all = fast_atof(value);
    } else if (same_str(key, "REFLECTIONS IN BIN    (WORKING SET)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().work_set_count = std::atoi(value);
    } else if (same_str(key, "BIN R VALUE           (WORKING SET)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_work = fast_atof(value);
    } else if (same_str(key, "BIN FREE R VALUE")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_free = fast_atof(value);
    } else if (same_str(key, "BIN FREE R VALUE TEST SET COUNT")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().rfree_set_count = std::atoi(value);
    } else if (same_str(key, "FROM WILSON PLOT           (A**2)")) {
      // TODO
      // exper.b_wilson = fast_atof(value);
    } else if (same_str(key, "MEAN B VALUE      (OVERALL, A**2)")) {
      ref_info.mean_b = fast_atof(value);
    } else if (same_str(key, "B11 (A**2)")) {
      ref_info.aniso_b.u11 = fast_atof(value);
    } else if (same_str(key, "B22 (A**2)")) {
      ref_info.aniso_b.u22 = fast_atof(value);
    } else if (same_str(key, "B33 (A**2)")) {
      ref_info.aniso_b.u33 = fast_atof(value);
    } else if (same_str(key, "B12 (A**2)")) {
      ref_info.aniso_b.u12 = fast_atof(value);
    } else if (same_str(key, "B13 (A**2)")) {
      ref_info.aniso_b.u13 = fast_atof(value);
    } else if (same_str(key, "B23 (A**2)")) {
      ref_info.aniso_b.u23 = fast_atof(value);
    } else if (same_str(key, "ESD FROM LUZZATI PLOT                    (A)")) {
      ref_info.luzzati_error = fast_atof(value);
    } else if (same_str(key, "DPI (BLOW EQ-10) BASED ON R VALUE        (A)")) {
      ref_info.dpi_blow_r = fast_atof(value);
    } else if (same_str(key, "DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A)")) {
      ref_info.dpi_blow_rfree = fast_atof(value);
    } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON R VALUE       (A)")) {
      ref_info.dpi_cruickshank_r = fast_atof(value);
    } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A)")) {
      ref_info.dpi_cruickshank_rfree = fast_atof(value);
    } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC")) {
      ref_info.cc_fo_fc_work = fast_atof(value);
    } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC FREE")) {
      ref_info.cc_fo_fc_free = fast_atof(value);
    } else if (same_str(key, "BOND LENGTHS")) {
      add_restraint_count_weight(ref_info, "t_bond_d", value);
    } else if (same_str(key, "BOND ANGLES")) {
      add_restraint_count_weight(ref_info, "t_angle_deg", value);
    } else if (same_str(key, "TORSION ANGLES")) {
      add_restraint_count_weight(ref_info, "t_dihedral_angle_d", value);
    } else if (same_str(key, "TRIGONAL CARBON PLANES")) {
      add_restraint_count_weight(ref_info, "t_trig_c_planes", value);
    } else if (same_str(key, "GENERAL PLANES")) {
      add_restraint_count_weight(ref_info, "t_gen_planes", value);
    } else if (same_str(key, "ISOTROPIC THERMAL FACTORS")) {
      add_restraint_count_weight(ref_info, "t_it", value);
    } else if (same_str(key, "BAD NON-BONDED CONTACTS")) {
      add_restraint_count_weight(ref_info, "t_nbd", value);
    } else if (same_str(key, "IMPROPER TORSIONS")) {
      add_restraint_count_weight(ref_info, "t_improper_torsion", value);
    } else if (same_str(key, "CHIRAL IMPROPER TORSION")) {
      add_restraint_count_weight(ref_info, "t_chiral_improper_torsion", value);
    } else if (same_str(key, "SUM OF OCCUPANCIES")) {
      add_restraint_count_weight(ref_info, "t_sum_occupancies", value);
    } else if (same_str(key, "UTILITY DISTANCES")) {
      add_restraint_count_weight(ref_info, "t_utility_distance", value);
    } else if (same_str(key, "UTILITY ANGLES")) {
      add_restraint_count_weight(ref_info, "t_utility_angle", value);
    } else if (same_str(key, "UTILITY TORSION")) {
      add_restraint_count_weight(ref_info, "t_utility_torsion", value);
    } else if (same_str(key, "IDEAL-DIST CONTACT TERM")) {
      add_restraint_count_weight(ref_info, "t_ideal_dist_contact", value);
    } else if (same_str(key, "BOND LENGTHS                       (A)")) {
      impl::find_or_add(ref_info.restr_stats, "t_bond_d").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "BOND ANGLES                  (DEGREES)")) {
      impl::find_or_add(ref_info.restr_stats, "t_angle_deg").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "PEPTIDE OMEGA TORSION ANGLES (DEGREES)")) {
      impl::find_or_add(ref_info.restr_stats, "t_omega_torsion").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "OTHER TORSION ANGLES         (DEGREES)")) {
      impl::find_or_add(ref_info.restr_stats, "t_other_torsion").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "TLS GROUP")) {
      ref_info.tls_groups.emplace_back();
      TlsGroup& tls_group = ref_info.tls_groups.back();
      tls_group.id = std::string(value, end);
      tls_group.num_id = (short) no_sign_atoi(tls_group.id.c_str());
    } else if (same_str(key, "SET") ||
               // "REMARK   3    SELECTION:"            -> TLS
               // "REMARK   3     SELECTION          :" -> NCS
               (same_str(key, "SELECTION") && colon == line + 23)) {
      if (!ref_info.tls_groups.empty()) {
        TlsGroup& group = ref_info.tls_groups.back();
        group.selections.emplace_back();
        group.selections.back().details = std::string(value, end);
        possibly_unfinished_remark3 = &group.selections.back().details;
      }
    } else if (same_str(key, "RESIDUE RANGE")) {
      if (!ref_info.tls_groups.empty() && end > colon+21) {
        TlsGroup& group = ref_info.tls_groups.back();
        group.selections.emplace_back();
        TlsGroup::Selection& sel = group.selections.back();
        sel.chain = read_string(colon+1, 5);
        if (sel.chain == read_string(colon+16, 5)) {
          try {
            sel.res_begin = SeqId(read_string(colon+6, 6));
            sel.res_end = SeqId(read_string(colon+21, 6));
          } catch (std::invalid_argument&) {
            group.selections.pop_back();
          }
        } else {  // unexpected -- TLS group should be in one chain
          group.selections.pop_back();
        }
      }
    } else if (same_str(key, "ORIGIN FOR THE GROUP (A)")) {
      std::vector<std::string> xyz = split_str_multi(std::string(value, end));
      if (ref_info.tls_groups.empty() || xyz.size() != 3)
        return;
      Position& origin = ref_info.tls_groups.back().origin;
      origin.x = fast_atof(xyz[0].c_str());
      origin.y = fast_atof(xyz[1].c_str());
      origin.z = fast_atof(xyz[2].c_str());
    } else if (is_tls_item(key)) {
      if (ref_info.tls_groups.empty())
        return;
      TlsGroup& tls = ref_info.tls_groups.back();
      std::vector<std::string> tokens = split_str_multi(key_start);
      for (size_t i = 0; i + 1 < tokens.size(); i += 2) {
        std::string& k = tokens[i];
        if (k.size() == 4 && k[3] == ':')
          k.resize(3);
        if (is_tls_item(k)) {
          int x = k[1] - '1';
          int y = k[2] - '1';
          double v = fast_atof(tokens[i+1].c_str());
          if (k[0] == 'S') {
            tls.S[x][y] = v;
          } else {
            SMat33<double>& tensor = k[0] == 'T' ? tls.T : tls.L;
            tensor.unchecked_ref(x, y) = v;
          }
        }
      }
    }
  } else {
    if (same_str(key, "DATA USED IN REFINEMENT.")) {
      meta.refinement.emplace_back();
      meta.refinement.back().id = std::to_string(meta.refinement.size());
    } else if (same_str(key, "FIT IN THE HIGHEST RESOLUTION BIN.")) {
      if (!meta.refinement.empty())
        meta.refinement.back().bins.emplace_back();
    }
  }
}

void read_remark_200_230_240(const char* line, Metadata& meta, std::string*& cryst_desc) {
  // multi-line continuation requires special handling
  if (cryst_desc) {
    if (line[10] == ' ' && line[11] == ' ') {
      const char* start = line + 11;
      cryst_desc->append(start, rtrim_cstr(start) - start);
      return;
    }
    cryst_desc = nullptr;
  }

  const char* key_start = skip_blank(line + 10);
  const char* colon = std::strchr(key_start, ':');
  const char* key_end = rtrim_cstr(key_start, colon);
  std::string key(key_start, key_end);
  if (colon) {
    const char* value = skip_blank(colon + 1);
    const char* end = rtrim_cstr(value);
    if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
      return;
    if (same_str(key, "INTENSITY-INTEGRATION SOFTWARE")) {
      add_software(meta, SoftwareItem::DataReduction, std::string(value, end));
    } else if (same_str(key, "DATA SCALING SOFTWARE")) {
      add_software(meta, SoftwareItem::DataScaling, std::string(value, end));
    } else if (same_str(key, "SOFTWARE USED")) {
      add_software(meta, SoftwareItem::Phasing, std::string(value, end));
    } else if (same_str(key, "METHOD USED TO DETERMINE THE STRUCTURE")) {
      meta.solved_by = std::string(value, end);
    } else if (same_str(key, "STARTING MODEL")) {
      meta.starting_model = std::string(value, end);
    } else if (!meta.experiments.empty()) {
      ExperimentInfo& exper = meta.experiments.back();
      DiffractionInfo& diffr = meta.crystals.back().diffractions[0];
      if (same_str(key, "EXPERIMENT TYPE")) {
        exper.method = std::string(value, end);
      } else if (same_str(key, "NUMBER OF CRYSTALS USED")) {
        exper.number_of_crystals = std::atoi(value);
      } else if (same_str(key, "PH")) {
        if (is_double(value))
          meta.crystals.back().ph = fast_atof(value);
        else
          meta.crystals.back().ph_range = std::string(value, end);
      } else if (same_str(key, "DATE OF DATA COLLECTION")) {
        diffr.collection_date = pdb_date_format_to_iso(std::string(value, end));
      } else if (same_str(key, "TEMPERATURE           (KELVIN)")) {
        diffr.temperature = fast_atof(value);
      } else if (same_str(key, "SYNCHROTRON              (Y/N)")) {
        if (*value == 'Y')
          diffr.source = "SYNCHROTRON";
      } else if (same_str(key, "RADIATION SOURCE")) {
        if (same_str(diffr.source, "SYNCHROTRON"))
          diffr.synchrotron = std::string(value, end);
        else
          diffr.source = std::string(value, end);
      } else if (same_str(key, "NEUTRON SOURCE")) {
        diffr.source = std::string(value, end);
      } else if (same_str(key, "BEAMLINE")) {
        diffr.beamline = std::string(value, end);
        if (!diffr.synchrotron.empty() && diffr.source_type.empty())
          diffr.source_type = diffr.synchrotron + " BEAMLINE " + diffr.beamline;
      } else if (same_str(key, "X-RAY GENERATOR MODEL")) {
        diffr.source_type = std::string(value, end);
      } else if (same_str(key, "MONOCHROMATIC OR LAUE    (M/L)")) {
        diffr.mono_or_laue = *value;
      } else if (same_str(key, "WAVELENGTH OR RANGE        (A)")) {
        diffr.wavelengths = std::string(value, end);
      } else if (same_str(key, "MONOCHROMATOR")) {
        diffr.monochromator = std::string(value, end);
      } else if (same_str(key, "OPTICS")) {
        diffr.optics = std::string(value, end);
      } else if (same_str(key, "DETECTOR TYPE")) {
        diffr.detector = std::string(value, end);
      } else if (same_str(key, "DETECTOR MANUFACTURER")) {
        diffr.detector_make = std::string(value, end);
      } else if (same_str(key, "NUMBER OF UNIQUE REFLECTIONS")) {
        exper.unique_reflections = std::atoi(value);
      } else if (same_str(key, "RESOLUTION RANGE HIGH      (A)")) {
        exper.reflections.resolution_high = fast_atof(value);
      } else if (same_str(key, "RESOLUTION RANGE LOW       (A)")) {
        exper.reflections.resolution_low = fast_atof(value);
      } else if (same_str(key, "COMPLETENESS FOR RANGE     (%)")) {
        exper.reflections.completeness = fast_atof(value);
      } else if (same_str(key, "DATA REDUNDANCY")) {
        exper.reflections.redundancy = fast_atof(value);
      } else if (same_str(key, "R MERGE                    (I)")) {
        exper.reflections.r_merge = fast_atof(value);
      } else if (same_str(key, "R SYM                      (I)")) {
        exper.reflections.r_sym = fast_atof(value);
      } else if (same_str(key, "<I/SIGMA(I)> FOR THE DATA SET")) {
        exper.reflections.mean_I_over_sigma = fast_atof(value);
      } else if (same_str(key, "REMARK")) {
        cryst_desc = &meta.crystals.back().description;
        *cryst_desc = std::string(value, end);
      } else if (!exper.shells.empty()) {
        if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)")) {
          exper.shells.back().resolution_high = fast_atof(value);
        } else if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE LOW  (A)")) {
          exper.shells.back().resolution_low = fast_atof(value);
        } else if (same_str(key, "COMPLETENESS FOR SHELL     (%)")) {
          exper.shells.back().completeness = fast_atof(value);
        } else if (same_str(key, "DATA REDUNDANCY IN SHELL")) {
          exper.shells.back().redundancy = fast_atof(value);
        } else if (same_str(key, "R MERGE FOR SHELL          (I)")) {
          exper.shells.back().r_merge = fast_atof(value);
        } else if (same_str(key, "R SYM FOR SHELL            (I)")) {
          exper.shells.back().r_sym = fast_atof(value);
        } else if (same_str(key, "<I/SIGMA(I)> FOR SHELL")) {
          exper.shells.back().mean_I_over_sigma = fast_atof(value);
        }
      }
    }
  } else {
    if (same_str(key, "EXPERIMENTAL DETAILS")) {
      meta.crystals.emplace_back();
      CrystalInfo& c = meta.crystals.back();
      c.id = std::to_string(meta.crystals.size());
      c.diffractions.emplace_back();
      c.diffractions[0].id = c.id;
      meta.experiments.emplace_back();
      meta.experiments.back().diffraction_ids.push_back(c.id);
      if (line[8] == '0' && line[9] == '0')
        c.diffractions[0].scattering_type = "x-ray";
      else if (line[8] == '3' && line[9] == '0')
        c.diffractions[0].scattering_type = "neutron";
      else if (line[8] == '4' && line[9] == '0')
        c.diffractions[0].scattering_type = "electron";
    }
    if (same_str(key, "IN THE HIGHEST RESOLUTION SHELL.")) {
      if (!meta.experiments.empty())
        meta.experiments.back().shells.emplace_back();
    }
  }
}

// Atom name and altloc are not provided in the SSBOND record.
// Usually it is SG (cysteine), but other disulfide bonds are also possible.
// If it's not SG, we pick the first sulfur atom in the residue.
const Residue* complete_ssbond_atom(AtomAddress& ad, const Model& mdl) {
  ad.atom_name = "SG";
  const_CRA cra = mdl.find_cra(ad);
  if (cra.residue && (!cra.atom || cra.atom->element != El::S))
    if (const Atom* a = cra.residue->find_by_element(El::S)) {
      ad.atom_name = a->name;
      ad.altloc = a->altloc;
    }
  return cra.residue;
}
void complete_ssbond(Connection& con, const Model& mdl, const UnitCell& cell) {
  const Residue* res1 = complete_ssbond_atom(con.partner1, mdl);
  const Residue* res2 = complete_ssbond_atom(con.partner2, mdl);
  if (res1 && res2 && (con.partner1.altloc != '\0' || con.partner2.altloc != '\0')) {
    // pick a pair of atoms in the shortest distance
    double min_dist_sq = INFINITY;
    for (const Atom& a1 : const_cast<Residue*>(res1)->get(con.partner1.atom_name))
      for (const Atom& a2 : const_cast<Residue*>(res2)->get(con.partner2.atom_name))
        if (a2.same_conformer(a1)) {
          double dist_sq = cell.find_nearest_image(a1.pos, a2.pos, con.asu).dist_sq;
          if (dist_sq < min_dist_sq) {
            con.partner1.altloc = a1.altloc;
            con.partner2.altloc = a2.altloc;
            min_dist_sq = dist_sq;
          }
        }
  }
}

Asu compare_link_symops(const std::string& record, short* reported_sym) {
  if (record.size() < 72)
    return Asu::Any;  // it could be interpreted as Same
  std::string s1 = read_string(&record[59], 6);
  std::string s2 = read_string(&record[66], 6);
  if (s1 == s2)
    return Asu::Same;
  size_t len1 = s1.length();
  size_t len2 = s2.length();
  if (len1 >= 4 && len1 < 6 && len2 >= 4 && len2 < 6) {
    // for 5 digits, we assume here that two digits are for sym_idx
    if (s1[0] == '1' && len1 == 4)  // symop1 is usually 1555
      reported_sym[0] = (short) read_int(s2.c_str(), len2 - 3);
    else
      reported_sym[0] = 99;
    for (size_t i = 1; i <= 3; ++i)
      reported_sym[i] = s2[len2 - 4 + i] - s1[len1 - 4 + i];
  }
  return Asu::Different;
}

void process_conn(Structure& st, const std::vector<std::string>& conn_records) {
  int disulf_count = 0;
  int covale_count = 0;
  int metalc_count = 0;
  for (const std::string& record : conn_records) {
    if (record[0] == 'S' || record[0] == 's') { // SSBOND
      if (record.length() < 32)
        continue;
      Connection c;
      c.name = "disulf" + std::to_string(++disulf_count);
      c.type = Connection::Disulf;
      const char* r = record.c_str();
      c.partner1.chain_name = read_string(r + 14, 2);
      c.partner1.res_id = read_res_id(r + 17, r + 11);
      c.partner2.chain_name = read_string(r + 28, 2);
      char res_id2[5] = {' ', ' ', ' ', ' ', ' '};
      std::memcpy(res_id2, r + 31, std::min((size_t)5, record.length() - 31));
      c.partner2.res_id = read_res_id(res_id2, r + 25);
      c.asu = compare_link_symops(record, c.reported_sym);
      if (record.length() > 73)
        c.reported_distance = read_double(r + 73, 5);
      complete_ssbond(c, st.first_model(), st.cell);
      st.connections.emplace_back(c);
    } else if (record[0] == 'L' || record[0] == 'l') { // LINK
      if (record.length() < 57)
        continue;
      Connection c;
      for (int i : {0, 1}) {
        const char* t = record.c_str() + 30 * i;
        AtomAddress& ad = (i == 0 ? c.partner1 : c.partner2);
        ad.chain_name = read_string(t + 20, 2);
        ad.res_id = read_res_id(t + 22, t + 17);
        ad.atom_name = read_string(t + 12, 4);
        ad.altloc = read_altloc(t[16]);
      }
      auto get_elem = [&](const char* name, const AtomAddress& ad) {
        if (element_from_padded_name_is_ambiguous(name)) {
          const_CRA cra = st.first_model().find_cra(ad);
          if (cra.atom)
            return cra.atom->element.elem;
        }
        return infer_element_from_padded_name(name);
      };
      // emulating names used in wwPDB mmCIFs (covaleN and metalcN)
      if (is_metal(get_elem(&record[12], c.partner1)) ||
          is_metal(get_elem(&record[42], c.partner2))) {
        c.name = "metalc" + std::to_string(++metalc_count);
        c.type = Connection::MetalC;
      } else {
        c.name = "covale" + std::to_string(++covale_count);
        c.type = Connection::Covale;
      }
      c.asu = compare_link_symops(record, c.reported_sym);
      if (record.length() > 73) {
        if (record[4] == 'R')
          c.link_id = read_string(&record[72], 8);
        else
          c.reported_distance = read_double(&record[73], 5);
      }
      st.connections.emplace_back(c);
    } else if (record[0] == 'C' || record[0] == 'c') { // CISPEP
      if (record.length() < 22)
        continue;
      const char* r = record.c_str();
      CisPep cispep;
      cispep.partner_c.chain_name = read_string(r + 14, 2);
      cispep.partner_c.res_id = read_res_id(r + 17, r + 11);
      cispep.partner_n.chain_name = read_string(r + 28, 2);
      cispep.partner_n.res_id = read_res_id(r + 31, r + 25);
      // In files with a single model in the PDB CISPEP modNum is 0,
      // but _struct_mon_prot_cis.pdbx_PDB_model_num is 1.
      cispep.model_num = st.models.size() == 1 ? st.models[0].num : read_int(r + 43, 3);
      cispep.reported_angle = read_double(r + 53, 6);
      st.cispeps.push_back(cispep);
    }
  }
}

// move initials after comma, as in mmCIF (A.-B.DOE -> DOE, A.-B.), see
// https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#AUTHOR
void change_author_name_format_to_mmcif(std::string& name) {
  // If the AUTHOR record has comma followed by space we get leading space here
  while (name[0] == ' ')
    name.erase(name.begin());
  size_t pos = 0;
  // Initials may have multiple letters (e.g. JU. or PON.)
  // but should not have space after dot.
  for (size_t i = 1; i < pos+4 && i+1 < name.size(); ++i)
    if (name[i] == '.' && name[i+1] != ' ')
      pos = i+1;
  if (pos > 0)
    name = name.substr(pos) + ", " + name.substr(0, pos);
}

// interprets subset of REMARKs from raw_remarks, filling in Metadata.
void read_metadata_from_remarks(Structure& st) {
  std::string* possibly_unfinished_remark3 = nullptr;
  std::string* cr_desc = nullptr;
  Transform matrix;
  for (const std::string& remark : st.raw_remarks) {
    if (remark.size() <= 11)
      continue;
    const char* line = remark.c_str();
    int num = read_int(line + 7, 3);
    switch (num) {
      case 2:
        if (st.resolution == 0.0 && std::strstr(line, "ANGSTROM"))
          st.resolution = read_double(line + 23, 7);
        break;
      case 3:
        read_remark3_line(line, st.meta, possibly_unfinished_remark3);
        break;
      case 200:
      case 230:
      case 240:
        read_remark_200_230_240(line, st.meta, cr_desc);
        break;
      case 300:
        if (!st.meta.remark_300_detail.empty()) {
          st.meta.remark_300_detail += '\n';
          st.meta.remark_300_detail += rtrim_str(remark.substr(11));
        } else if (remark.compare(11, 7, "REMARK:") == 0) {
          st.meta.remark_300_detail = trim_str(remark.substr(18));
        }
        break;
      case 350: {
        const char* colon = std::strchr(line+11, ':');
        if (colon == line+22 && starts_with(line+11, "BIOMOLECULE")) {
          st.assemblies.emplace_back(read_string(line+23, 20));
          continue;
        }
        if (st.assemblies.empty())
          continue;
        Assembly& assembly = st.assemblies.back();
        auto r350_key = [&](int cpos, const char* text) {
          return colon == line + cpos && starts_with(line+11, text);
        };
        if (starts_with(line+11, "  BIOMT")) {
          if (read_matrix(matrix, line+13, remark.size()-13) == 3)
            if (!assembly.generators.empty()) {
              auto& opers = assembly.generators.back().operators;
              opers.emplace_back();
              opers.back().name = read_string(line+20, 3);
              opers.back().transform = matrix;
              matrix.set_identity();
            }
        } else if (r350_key(44, "AUTHOR DETERMINED")) {
          assembly.author_determined = true;
          assembly.oligomeric_details = read_string(line+45, 35);
        } else if (r350_key(51, "SOFTWARE DETERMINED")) {
          assembly.software_determined = true;
          assembly.oligomeric_details = read_string(line+52, 28);
        } else if (r350_key(24, "SOFTWARE USED")) {
          assembly.software_name = read_string(line+25, 55);
        } else if (r350_key(36, "TOTAL BURIED SURFACE AREA")) {
          assembly.absa = read_double(line+37, 12);
        } else if (r350_key(38, "SURFACE AREA OF THE COMPLEX")) {
          assembly.ssa = read_double(line+39, 12);
        } else if (r350_key(40, "CHANGE IN SOLVENT FREE ENERGY")) {
          assembly.more = read_double(line+41, 12);
        } else if (r350_key(40, "APPLY THE FOLLOWING TO CHAINS") ||
                   r350_key(40, "                   AND CHAINS")) {
          if (line[11] == 'A') // first line - APPLY ...
            assembly.generators.emplace_back();
          else if (assembly.generators.empty())
            continue;
          split_str_into_multi(read_string(line+41, 39), ", ",
                               assembly.generators.back().chains);
        }
      }
    }
    // if REMARK 2 was missing, try resolution from REMARK 3
    if (st.resolution == 0.0) {
      for (const RefinementInfo& ref_info : st.meta.refinement)
        if (!std::isnan(ref_info.resolution_high) && ref_info.resolution_high != 0.) {
          st.resolution = ref_info.resolution_high;
          break;
        }
    }
  }
}

} // anonymous namespace

Structure read_pdb_from_stream(AnyStream& line_reader, const std::string& source,
                               PdbReadOptions options) {
  if (options.max_line_length <= 0 || options.max_line_length > 120)
    options.max_line_length = 120;
  Structure st;
  st.input_format = CoorFormat::Pdb;
  st.name = path_basename(source, {".gz", ".pdb"});
  Transform matrix;
  std::vector<std::string> conn_records;
  std::unordered_map<ResidueId, int> resmap;
  Model *model = nullptr;
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  char line[122] = {0};
  int line_num = 0;
  bool after_ter = false;
  auto wrong = [&line_num](const std::string& msg) {
    fail("Problem in line ", std::to_string(line_num), ": ", msg);
  };
  while (size_t len = line_reader.copy_line(line, options.max_line_length+1)) {
    ++line_num;
    if (is_record_type4(line, "ATOM") || is_record_type4(line, "HETATM")) {
      if (len < 55)
        wrong("The line is too short to be correct:\n" + std::string(line));
      std::string chain_name = read_string(line+20, 2);
      ResidueId rid = read_res_id(line+22, line+17);

      if (!chain || chain_name != chain->name) {
        if (!model) {
          // A single model usually doesn't have the MODEL record. Also,
          // MD trajectories may have frames separated by ENDMDL without MODEL.
          int num = (int) st.models.size() + 1;
          if (st.find_model(num))
            wrong("ATOM/HETATM between models");
          st.models.emplace_back(num);
          model = &st.models.back();
        }
        const Chain* prev_part = model->find_chain(chain_name);
        after_ter = prev_part &&
                    prev_part->residues[0].entity_type == EntityType::Polymer;
        model->chains.emplace_back(chain_name);
        chain = &model->chains.back();
        resmap.clear();
        resi = nullptr;
      }
      // Non-standard but widely used 4-character segment identifier.
      // Left-justified, and may include a space in the middle.
      // The segment may be a portion of a chain or a complete chain.
      if (len > 72)
        rid.segment = read_string(line+72, 4);
      if (!resi || !resi->matches(rid)) {
        auto it = resmap.find(rid);
        // In normal PDB files it is fast enough to use
        // resi = chain->find_residue(rid);
        // but in pseudo-PDB files (such as MD files where millions
        // of residues are in the same "chain") it is too slow.
        if (it == resmap.end()) {
          resmap.emplace(rid, (int) chain->residues.size());
          chain->residues.emplace_back(rid);
          resi = &chain->residues.back();

          resi->het_flag = line[0] & ~0x20;
          if (after_ter)
            resi->entity_type = resi->is_water() ? EntityType::Water
                                                 : EntityType::NonPolymer;
        } else {
          resi = &chain->residues[it->second];
        }
      }

      Atom atom;
      atom.serial = read_serial(line+6);
      atom.name = read_string(line+12, 4);
      atom.altloc = read_altloc(line[16]);
      atom.pos.x = read_double(line+30, 8);
      atom.pos.y = read_double(line+38, 8);
      atom.pos.z = read_double(line+46, 8);
      if (len > 58)
        atom.occ = (float) read_double(line+54, 6);
      if (len > 64)
        atom.b_iso = (float) read_double(line+60, 6);
      if (len > 76 && (std::isalpha(line[76]) || std::isalpha(line[77])))
        atom.element = Element(line + 76);
      else
        atom.element = infer_element_from_padded_name(line+12);
      atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
      resi->atoms.emplace_back(atom);

    } else if (is_record_type4(line, "ANISOU")) {
      if (!model || !chain || !resi || resi->atoms.empty())
        wrong("ANISOU record not directly after ATOM/HETATM.");
      // We assume that ANISOU refers to the last atom.
      // Can it not be the case?
      Atom &atom = resi->atoms.back();
      if (atom.aniso.u11 != 0.)
        wrong("Duplicated ANISOU record or not directly after ATOM/HETATM.");
      atom.aniso.u11 = read_int(line+28, 7) * 1e-4f;
      atom.aniso.u22 = read_int(line+35, 7) * 1e-4f;
      atom.aniso.u33 = read_int(line+42, 7) * 1e-4f;
      atom.aniso.u12 = read_int(line+49, 7) * 1e-4f;
      atom.aniso.u13 = read_int(line+56, 7) * 1e-4f;
      atom.aniso.u23 = read_int(line+63, 7) * 1e-4f;

    } else if (is_record_type4(line, "REMARK")) {
      if (line[len-1] == '\n')
        --len;
      if (line[len-1] == '\r')
        --len;
      st.raw_remarks.emplace_back(line, line+len);

    } else if (is_record_type4(line, "CONECT")) {
      int serial = read_serial(line+6);
      if (len >= 11 && serial != 0) {
        std::vector<int>& bonded_atoms = st.conect_map[serial];
        int limit = std::min(27, (int)len - 1);
        for (int offset = 11; offset <= limit; offset += 5) {
          int n = read_serial(line+offset);
          if (n != 0)
            bonded_atoms.push_back(n);
        }
      }

    } else if (is_record_type4(line, "SEQRES")) {
      std::string chain_name = read_string(line+10, 2);
      Entity& ent = impl::find_or_add(st.entities, chain_name);
      ent.entity_type = EntityType::Polymer;
      for (int i = 19; i < 68; i += 4) {
        std::string res_name = read_string(line+i, 3);
        if (!res_name.empty())
          ent.full_sequence.emplace_back(res_name);
      }

    } else if (is_record_type4(line, "HELIX")) {
      if (len < 40)
        continue;
      Helix helix;
      helix.start.chain_name = read_string(line+18, 2);
      helix.start.res_id = read_res_id(line+21, line+15);
      helix.end.chain_name = read_string(line+30, 2);
      helix.end.res_id = read_res_id(line+33, line+27);
      helix.set_helix_class_as_int(read_int(line+38, 2));
      if (len > 72)
        helix.length = read_int(line+72, 5);
      st.helices.emplace_back(helix);

    } else if (is_record_type4(line, "SHEET")) {
      if (len < 40)
        continue;
      std::string sheet_id = read_string(line+11, 3);
      Sheet& sheet = impl::find_or_add(st.sheets, sheet_id);
      sheet.strands.emplace_back();
      Sheet::Strand& strand = sheet.strands.back();
      strand.start.chain_name = read_string(line+20, 2);
      strand.start.res_id = read_res_id(line+22, line+17);
      strand.end.chain_name = read_string(line+31, 2);
      strand.end.res_id = read_res_id(line+33, line+28);
      strand.sense = read_int(line+38, 2);
      if (len > 67) {
        // the SHEET record has no altloc for atoms of hydrogen bond
        strand.hbond_atom2.atom_name = read_string(line+41, 4);
        strand.hbond_atom2.chain_name = read_string(line+48, 2);
        strand.hbond_atom2.res_id = read_res_id(line+50, line+45);
        strand.hbond_atom1.atom_name = read_string(line+56, 4);
        strand.hbond_atom1.chain_name = read_string(line+63, 2);
        strand.hbond_atom1.res_id = read_res_id(line+65, line+60);
      }

    } else if (is_record_type4(line, "SSBOND") ||
               is_record_type4(line, "LINK") ||
               is_record_type4(line, "CISPEP")) {
      conn_records.emplace_back(line);

    } else if (is_record_type3(line, "TER")) { // finishes polymer chains
      if (!chain || st.ter_status == 'e')
        continue;
      st.ter_status = 'y';
      if (options.split_chain_on_ter) {
        chain = nullptr;
        // split_chain_on_ter is used for AMBER files that can have TER records
        // in various places. So in such case TER doesn't imply entity_type.
        continue;
      }
      // If we have 2+ TER records in one chain, they are used in non-standard
      // way and should be better ignored (in all the chains).
      if (after_ter) {
        st.ter_status = 'e';  // all entity_types will be later set to Unknown
        continue;
      }
      for (Residue& res : chain->residues) {
        res.entity_type = EntityType::Polymer;
        // Sanity check: water should not be marked as a polymer.
        if GEMMI_UNLIKELY(res.is_water())
          st.ter_status = 'e';  // all entity_types will be later set to Unknown
      }
      after_ter = true;

    } else if (is_record_type4(line, "MODRES")) {
      ModRes modres;
      modres.chain_name = read_string(line + 15, 2);
      modres.res_id = read_res_id(line + 18, line + 12);
      modres.parent_comp_id = read_string(line + 24, 3);
      if (len >= 30)
        // this field is named comment in PDB spec, but details in mmCIF
        modres.details = read_string(line + 29, 41);
      // Refmac's extension: 73-80 mod_id
      // Check for spaces to make sure it's not an overflowed comment
      if (len >= 73 && line[70] == ' ' && line[71] == ' ')
        modres.mod_id = read_string(line + 72, 8);
      st.mod_residues.push_back(modres);

    } else if (is_record_type4(line, "HETNAM")) {
      if (len > 71 && line[70] == ' ') {
        std::string full_code = read_string(line + 71, 8);
        if (!full_code.empty())
          st.shortened_ccd_codes.emplace_back(full_code, read_string(line + 11, 3));
      }

    } else if (is_record_type4(line, "DBREF")) { // DBREF or DBREF1 or DBREF2
      std::string chain_name = read_string(line+11, 2);
      Entity& ent = impl::find_or_add(st.entities, chain_name);
      ent.entity_type = EntityType::Polymer;
      if (line[5] == ' ' || line[5] == '1')
        ent.dbrefs.emplace_back();
      else if (ent.dbrefs.empty()) // DBREF2 without DBREF1?
        break;
      Entity::DbRef& dbref = ent.dbrefs.back();
      if (line[5] == ' ' || line[5] == '1') {
        dbref.seq_begin = read_seq_id(line+14);
        dbref.seq_end = read_seq_id(line+20);
        dbref.db_name = read_string(line+26, 6);
        if (line[5] == ' ') {
          dbref.accession_code = read_string(line+33, 8);
          dbref.id_code = read_string(line+42, 12);
          dbref.db_begin.num = read_int(line+55, 5);
          dbref.db_begin.icode = line[60];
          dbref.db_end.num = read_int(line+62, 5);
          dbref.db_end.icode = line[67];
        } else {  // line[5] == '1'
          dbref.id_code = read_string(line+47, 20);
        }
      } else if (line[5] == '2') {
        dbref.accession_code = read_string(line+18, 22);
        dbref.db_begin.num = read_int(line+45, 10);
        dbref.db_end.num = read_int(line+57, 10);
      }

    } else if (is_record_type4(line, "HEADER")) {
      if (len > 50)
        st.info["_struct_keywords.pdbx_keywords"] = rtrim_str(std::string(line+10, 40));
      if (len > 59) { // date in PDB has format 28-MAR-07
        std::string date = pdb_date_format_to_iso(std::string(line+50, 9));
        if (!date.empty())
          st.info["_pdbx_database_status.recvd_initial_deposition_date"] = date;
      }
      if (len > 66) {
        std::string entry_id = rtrim_str(std::string(line+62, 4));
        if (!entry_id.empty())
          st.info["_entry.id"] = entry_id;
      }

    } else if (is_record_type4(line, "TITLE")) {
      if (len > 10)
        st.info["_struct.title"] += rtrim_str(std::string(line+10, len-10-1));

    } else if (is_record_type4(line, "KEYWDS")) {
      if (len > 10)
        st.info["_struct_keywords.text"] += rtrim_str(std::string(line+10, len-10-1));

    } else if (is_record_type4(line, "EXPDTA")) {
      if (len > 10)
        st.info["_exptl.method"] += trim_str(std::string(line+10, len-10-1));

    } else if (is_record_type4(line, "AUTHOR") && len > 10) {
      std::string last;
      if (!st.meta.authors.empty()) {
        last = st.meta.authors.back();
        st.meta.authors.pop_back();
      }
      size_t prev_size = st.meta.authors.size();
      const char* start = skip_blank(line+10);
      const char* end = rtrim_cstr(start, line+len);
      split_str_into(std::string(start, end), ',', st.meta.authors);
      if (!last.empty() && st.meta.authors.size() > prev_size) {
        // the spaces were trimmed, we may need a space between words
        if (last.back() != '-' && last.back() != '.')
          last += ' ';
        st.meta.authors[prev_size].insert(0, last);
      }

    } else if (is_record_type4(line, "SCALEn")) {
      if (read_matrix(matrix, line, len) == 3) {
        st.cell.set_matrices_from_fract(matrix);
        matrix.set_identity();
      }

    } else if (is_record_type4(line, "ORIGX")) {
      st.has_origx = true;
      read_matrix(st.origx, line, len);

    } else if (is_record_type4(line, "CRYST1")) {
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

    } else if (is_record_type4(line, "MTRIXn")) {
      if (read_matrix(matrix, line, len) == 3) {
        std::string id = read_string(line+7, 3);
        if (matrix.is_identity()) {
          // store only ID that will be used when writing to file
          st.info["_struct_ncs_oper.id"] = id;
        } else {
          bool given = len > 59 && line[59] == '1';
          st.ncs.push_back({id, given, matrix});
          matrix.set_identity();
        }
      }
    } else if (is_record_type4(line, "MODEL")) {
      if (model && chain)
        wrong("MODEL without ENDMDL?");
      int num = read_int(line+6, 8);
      model = &st.find_or_add_model(num);
      if (!model->chains.empty())
        wrong("duplicate MODEL number: " + std::to_string(num));
      chain = nullptr;

    } else if (is_record_type4(line, "ENDMDL")) {
      model = nullptr;
      chain = nullptr;

    } else if (is_record_type3(line, "END")) {
      break;
    } else if (is_record_type4(line, "data")) {
      if (line[4] == '_' && !model)
        fail("Incorrect file format (perhaps it is cif not pdb?): " + source);
    } else if (is_record_type4(line, "{\"da")) {
      if (ialpha3_id(line+4) == ialpha3_id("ta_") && !model)
      fail("Incorrect file format (perhaps it is mmJSON not pdb?): " + source);
    }
  }
  // If we read a PDB header (they can be downloaded from RSCB) we have no
  // models. User's code may not expect this. Usually, empty model will be
  // handled more gracefully than no models.
  if (st.models.empty())
    st.models.emplace_back(1);

  if (st.ter_status == 'e')
    remove_entity_types(st);

  // Here we assign Residue::subchain, but only for chains with all
  // Residue::entity_type assigned, i.e. for chains with TER.
  assign_subchains(st, /*force=*/false, /*fail_if_unknown=*/false);

  for (Chain& ch : st.models[0].chains)
    if (Entity* entity = st.get_entity(ch.name))
      if (auto polymer = ch.get_polymer())
        entity->subchains.emplace_back(polymer.subchain_id());

  st.setup_cell_images();

  process_conn(st, conn_records);

  for (std::string& name : st.meta.authors)
    change_author_name_format_to_mmcif(name);

  if (!options.skip_remarks)
    read_metadata_from_remarks(st);

  restore_full_ccd_codes(st);
  return st;
}

std::vector<Op> read_remark_290(const std::vector<std::string>& raw_remarks) {
  std::vector<Op> ops;
  // we only check triplet notation:
  // REMARK 290     NNNMMM   OPERATOR
  // REMARK 290       1555   X,Y,Z
  for (const std::string& remark : raw_remarks)
    if (remark.size() > 25 && std::memcmp(&remark[7], "290", 3) == 0 &&
        std::memcmp(&remark[10], "     ", 5) == 0 &&
        std::memcmp(&remark[18], "555   ", 6) == 0) {
      if (read_int(remark.c_str() + 15, 3) != (int)ops.size() + 1)
        fail("Symmetry operators not in order?: " + remark);
      Op op = parse_triplet(read_string(remark.c_str() + 24, 56));
      ops.push_back(op);
    }
  return ops;
}

} // namespace gemmi
