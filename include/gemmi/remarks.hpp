// Copyright 2019 Global Phasing Ltd.
//
// Reading some of the REMARKs from the PDB format.

#ifndef GEMMI_REMARKS_HPP_
#define GEMMI_REMARKS_HPP_

#include "metadata.hpp"
#include "pdb.hpp"  // for pdb_date_format_to_iso, read_double, ...

namespace gemmi {

namespace pdb_impl {

inline bool is_double(const char* p) {
  while (std::isspace(*p))
    ++p;
  if (*p == '-' || *p == '+')
    ++p;
  while (*p >= '0' && *p <= '9')
    ++p;
  if (*p == '.') {
    ++p;
    while (*p >= '0' && *p <= '9')
      ++p;
  }
  while (std::isspace(*p))
    ++p;
  return *p == '\0';
}

inline bool is_tls_item(const std::string& key) {
  return key.size() == 3 &&
    (key[0] == 'T' || key[0] == 'L' || key[0] == 'S') &&
    (key[1] == '1' || key[1] == '2' || key[1] == '3') &&
    (key[2] == '1' || key[2] == '2' || key[2] == '3');
}


inline void read_remark3_line(const char* line, Metadata& meta) {
  // Based on:
  // www.wwpdb.org/documentation/file-format-content/format23/remark3.html
  // and analysis of PDB files.
  // In special cases, such as joint X-ray and neutron refinement 5MOO,
  // PDB file can have two REMARK 3 blocks.
  // Generally, after "REMARK   3" we have either a header-like sentance
  // or a key:value pair with a colon, or a continuation of text from the
  // previous line.
  const char* key_start = skip_blank(line + 10);
  const char* colon = std::strchr(key_start, ':');
  const char* key_end = rtrim_cstr(key_start, colon);
  std::string key(key_start, key_end);
  if (colon) {
    const char* value = skip_blank(colon + 1);
    const char* end = rtrim_cstr(value);
    if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
      return;
    if (same_str(key, "PROGRAM"))
      SoftwareItem& item = meta.add_software(SoftwareItem::Refinement,
                                             std::string(value, end));
    if (meta.refinement.empty())
      return;
    RefinementInfo& ref_info = meta.refinement.back();
    if (same_str(key, "RESOLUTION RANGE HIGH (ANGSTROMS)")) {
      ref_info.resolution_high = read_double(value);
    } else if (same_str(key, "RESOLUTION RANGE LOW  (ANGSTROMS)")) {
      ref_info.resolution_low = read_double(value);
    } else if (same_str(key, "COMPLETENESS FOR RANGE        (%)")) {
      ref_info.completeness = read_double(value);
    } else if (same_str(key, "NUMBER OF REFLECTIONS")) {
      ref_info.reflection_count = std::atoi(value);
    } else if (same_str(key, "CROSS-VALIDATION METHOD")) {
      ref_info.cross_validation_method = std::string(value, end);
    } else if (same_str(key, "FREE R VALUE TEST SET SELECTION")) {
      ref_info.rfree_selection_method = std::string(value, end);
    } else if (same_str(key, "R VALUE     (WORKING + TEST SET)")) {
      ref_info.r_all = read_double(value);
    } else if (same_str(key, "R VALUE            (WORKING SET)")) {
      ref_info.r_work = read_double(value);
    } else if (same_str(key, "FREE R VALUE")) {
      ref_info.r_free = read_double(value);
    } else if (same_str(key, "FREE R VALUE TEST SET COUNT")) {
      ref_info.rfree_set_count = atoi(value);
    } else if (same_str(key, "TOTAL NUMBER OF BINS USED")) {
      ref_info.bin_count = std::atoi(value);
    } else if (same_str(key, "BIN RESOLUTION RANGE HIGH       (A)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().resolution_high = read_double(value);
    } else if (same_str(key, "BIN RESOLUTION RANGE LOW        (A)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().resolution_low = read_double(value);
    } else if (same_str(key, "BIN COMPLETENESS (WORKING+TEST) (%)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().completeness = read_double(value);
    } else if (same_str(key, "REFLECTIONS IN BIN   (WORKING+TEST)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().reflection_count = std::atoi(value);
    } else if (same_str(key, "BIN R VALUE          (WORKING+TEST)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_all = read_double(value);
    } else if (same_str(key, "BIN R VALUE           (WORKING SET)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_work = read_double(value);
    } else if (same_str(key, "BIN FREE R VALUE")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_free = read_double(value);
    } else if (same_str(key, "BIN FREE R VALUE TEST SET COUNT")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().rfree_set_count = std::atoi(value);
    } else if (same_str(key, "FROM WILSON PLOT           (A**2)")) {
      // TODO
      // exper.b_wilson = read_double(value);
    } else if (same_str(key, "MEAN B VALUE      (OVERALL, A**2)")) {
      ref_info.mean_b = read_double(value);
    } else if (same_str(key, "B11 (A**2)")) {
      ref_info.aniso_b[0][0] = read_double(value);
    } else if (same_str(key, "B22 (A**2)")) {
      ref_info.aniso_b[1][1] = read_double(value);
    } else if (same_str(key, "B33 (A**2)")) {
      ref_info.aniso_b[2][2] = read_double(value);
    } else if (same_str(key, "B12 (A**2)")) {
      ref_info.aniso_b[0][1] = read_double(value);
    } else if (same_str(key, "B13 (A**2)")) {
      ref_info.aniso_b[0][2] = read_double(value);
    } else if (same_str(key, "B23 (A**2)")) {
      ref_info.aniso_b[1][2] = read_double(value);
    } else if (same_str(key, "ESD FROM LUZZATI PLOT                    (A)")) {
      ref_info.luzzati_error = read_double(value);
    } else if (same_str(key, "DPI (BLOW EQ-10) BASED ON R VALUE        (A)")) {
      ref_info.dpi_blow_r = read_double(value);
    } else if (same_str(key, "DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A)")) {
      ref_info.dpi_blow_rfree = read_double(value);
    } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON R VALUE       (A)")) {
      ref_info.dpi_cruickshank_r = read_double(value);
    } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A)")) {
      ref_info.dpi_cruickshank_rfree = read_double(value);
    } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC")) {
      ref_info.cc_fo_fc = read_double(value);
    } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC FREE")) {
      ref_info.cc_fo_fc_free = read_double(value);
    } else if (same_str(key, "TLS GROUP")) {
      ref_info.tls_groups.emplace_back();
      ref_info.tls_groups.back().id = std::string(value, end);
    } else if (same_str(key, "SET")) {
      if (!ref_info.tls_groups.empty())
        ref_info.tls_groups.back().selection = std::string(value, end);
    } else if (same_str(key, "ORIGIN FOR THE GROUP (A)")) {
      std::vector<std::string> xyz = split_str_multi(std::string(value, end));
      if (ref_info.tls_groups.empty() || xyz.size() != 3)
        return;
      Position& origin = ref_info.tls_groups.back().origin;
      origin.x = read_double(xyz[0].c_str());
      origin.y = read_double(xyz[1].c_str());
      origin.z = read_double(xyz[2].c_str());
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
          Mat33& m = k[0] == 'T' ? tls.T : k[0] == 'L' ? tls.L : tls.S;
          int x = k[1] - '1';
          int y = k[2] - '1';
          m[x][y] = m[y][x] = read_double(tokens[i+1].c_str());
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

inline void read_remark_200_230_240(const char* line, Metadata& meta) {
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
      meta.add_software(SoftwareItem::DataReduction, std::string(value, end));
    } else if (same_str(key, "DATA SCALING SOFTWARE")) {
      meta.add_software(SoftwareItem::DataScaling, std::string(value, end));
    } else if (same_str(key, "SOFTWARE USED")) {
      meta.add_software(SoftwareItem::Phasing, std::string(value, end));
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
          meta.crystals.back().ph = read_double(value);
        else
          meta.crystals.back().ph_range = std::string(value, end);
      } else if (same_str(key, "DATE OF DATA COLLECTION")) {
        diffr.collection_date = pdb_date_format_to_iso(std::string(value, end));
      } else if (same_str(key, "TEMPERATURE           (KELVIN)")) {
        diffr.temperature = read_double(value);
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
        exper.reflections.resolution_high = read_double(value);
      } else if (same_str(key, "RESOLUTION RANGE LOW       (A)")) {
        exper.reflections.resolution_low = read_double(value);
      } else if (same_str(key, "COMPLETENESS FOR RANGE     (%)")) {
        exper.reflections.completeness = read_double(value);
      } else if (same_str(key, "DATA REDUNDANCY")) {
        exper.reflections.redundancy = read_double(value);
      } else if (same_str(key, "R MERGE                    (I)")) {
        exper.reflections.r_merge = read_double(value);
      } else if (same_str(key, "R SYM                      (I)")) {
        exper.reflections.r_sym = read_double(value);
      } else if (same_str(key, "<I/SIGMA(I)> FOR THE DATA SET")) {
        exper.reflections.mean_I_over_sigma = read_double(value);
      } else if (!exper.shells.empty()) {
        if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)")) {
          exper.shells.back().resolution_high = read_double(value);
        } else if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE LOW  (A)")) {
          exper.shells.back().resolution_low = read_double(value);
        } else if (same_str(key, "COMPLETENESS FOR SHELL     (%)")) {
          exper.shells.back().completeness = read_double(value);
        } else if (same_str(key, "DATA REDUNDANCY IN SHELL")) {
          exper.shells.back().redundancy = read_double(value);
        } else if (same_str(key, "R MERGE FOR SHELL          (I)")) {
          exper.shells.back().r_merge = read_double(value);
        } else if (same_str(key, "R SYM FOR SHELL            (I)")) {
          exper.shells.back().r_sym = read_double(value);
        } else if (same_str(key, "<I/SIGMA(I)> FOR SHELL")) {
          exper.shells.back().mean_I_over_sigma = read_double(value);
        }
      }
    }
  } else {
    if (same_str(key, "EXPERIMENTAL DETAILS")) {
      meta.experiments.emplace_back();
      meta.crystals.emplace_back();
      CrystalInfo& c = meta.crystals.back();
      c.id = std::to_string(meta.crystals.size());
      c.diffractions.emplace_back();
      c.diffractions[0].id = c.id;
      c.diffractions[0].crystal_id = c.id;
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

} // namespace pdb_impl

void read_metadata_from_remarks(Structure& st) {
  for (const std::string& remark : st.raw_remarks)
    if (remark.size() > 11) {
      switch (pdb_impl::read_int(remark.c_str() + 7, 3)) {
        case 3:
          pdb_impl::read_remark3_line(remark.c_str(), st.meta);
          break;
        case 200:
        case 230:
        case 240:
          pdb_impl::read_remark_200_230_240(remark.c_str(), st.meta);
          break;
      }
    }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et

