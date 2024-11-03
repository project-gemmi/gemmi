// Copyright 2019 Global Phasing Ltd.
//
// Transform MTZ or SF-mmCIF map coefficients to CCP4 map.

#include <stdio.h>
#include <unordered_map>
#include <gemmi/ccp4.hpp>      // for Ccp4
#include <gemmi/calculate.hpp> // for calculate_fractional_box
#include <gemmi/select.hpp>    // for Selection
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/mtz.hpp>       // for Mtz
#include "mapcoef.h"

#define GEMMI_PROG sf2map
#include "options.h"

namespace {

enum OptionIndex { Normalize=AfterMapOptions, MapMask, Margin, Select, Check };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE MAP_FILE\n"
    "  " EXE_NAME " --check INPUT_FILE\n\n"
    "INPUT_FILE must be an MTZ or mmCIF file with map coefficients.\n\n"
    "By default, the program searches for 2mFo-DFc map coefficients in:\n"
    "  - MTZ columns FWT/PHWT or 2FOFCWT/PH2FOFCWT,\n"
    "  - mmCIF tags _refln.pdbx_FWT/pdbx_PHWT.\n"
    "If option \"-d\" is given, mFo-DFc map coefficients are searched in:\n"
    "  - MTZ columns DELFWT/PHDELWT or FOFCWT/PHFOFCWT,\n"
    "  - mmCIF tags _refln.pdbx_DELFWT/pdbx_DELPHWT.\n\n"
    "Option --check can take column name mapping as arg, for example:\n"
    "  --check=dfc:FCALC,fosc:FOBS  (typically, use --check without an arg).\n"
    "Possible keys are fwt, phwt, delfwt, phdelwt, dfc, fosc, fom, free.\n"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  MapUsage[Diff],
  MapUsage[Section],
  MapUsage[FLabel],
  MapUsage[PhLabel],
  MapUsage[WeightLabel],
  MapUsage[GridDims],
  MapUsage[ExactDims],
  MapUsage[Sample],
  MapUsage[AxesZyx],
  MapUsage[GridQuery],
  MapUsage[TimingFft],
  { Normalize, 0, "", "normalize", Arg::None,
    "  --normalize  \tScale the map to standard deviation 1 and mean 0." },
  { MapMask, 0, "", "mapmask", Arg::Required,
    "  --mapmask=FILE  \tOutput only map covering the structure from FILE,"
    " similarly to CCP4 MAPMASK with XYZIN." },
  { Margin, 0, "", "margin", Arg::Float,
    "  --margin=N  \t(w/ --mapmask) Border in Angstroms (default: 5)." },
  { Select, 0, "", "select", Arg::Required,
    "  --select=SEL  \t(w/ --mapmask) Atom selection for mask, MMDB syntax." },
  { Check, 0, "", "check", Arg::Optional,
    "  --check[=cols]  \tAnalyze map coefficient columns in MTZ file." },
  { 0, 0, 0, 0, 0, 0 }
};


void transform_sf_to_map(OptParser& p) {
  const char* input_path = p.nonOption(0);
  const char* map_path = p.options[GridQuery] ? nullptr : p.nonOption(1);
  gemmi::Ccp4<float> ccp4;
  ccp4.grid = read_sf_and_fft_to_map(input_path, p.options,
                                     p.options[Verbose] ? stderr : nullptr);
  if (p.options[Verbose])
    fprintf(stderr, "Writing %s ...\n", map_path);
  ccp4.update_ccp4_header(2);
  if (p.options[Normalize]) {
    double mult = 1.0 / ccp4.hstats.rms;
    for (float& x : ccp4.grid.data)
      x = float((x - ccp4.hstats.dmean) * mult);
    ccp4.update_ccp4_header(2);
  }
  if (p.options[MapMask]) {
    double margin = 5;
    if (p.options[Margin])
      margin = std::atof(p.options[Margin].arg);
    gemmi::Structure st = gemmi::read_structure_gz(p.options[MapMask].arg);
    st.cell = ccp4.grid.unit_cell; // needed in case the unit cells differ
    gemmi::Box<gemmi::Fractional> box;
    if (p.options[Select]) {
      gemmi::Selection sel(p.options[Select].arg);
      for (gemmi::Model& model : sel.models(st))
        for (gemmi::Chain& chain : sel.chains(model))
          for (gemmi::Residue& res : sel.residues(chain))
            for (gemmi::Atom& atom : sel.atoms(res))
              box.extend(st.cell.fractionalize(atom.pos));
      box.add_margins({margin * st.cell.ar, margin * st.cell.br, margin * st.cell.cr});
    } else {
      box = gemmi::calculate_fractional_box(st, margin);
    }
    ccp4.set_extent(box);
  }
  ccp4.write_ccp4_map(map_path);
}

struct ReflData {
  gemmi::Miller hkl;
  int free_flag;
  float fm;     // amplitude from a normal, 2-1 map (FWT, 2FOFCWT or similar)
  float fm_phi; // phase for fm, in degrees
  float fd;     // amplitude of a difference map (DELFWT, FOFCWT or similar)
  float fd_phi; // phase for fd
  float fc;     // D*Fc w/ bulk solvent, probably
  float phi;    // calculated phase (it's used also for Fo, which has no phase)
  float fo;     // scaled Fo, probably
  float m;      // figure of merit (FOM)
};

struct DataForCheck {
  const gemmi::SpaceGroup* sg;
  bool has_free_flags;
  bool has_fo;
  std::vector<ReflData> refl;
  std::array<std::string, 9> labels;  // corresponds to ReflData items, except hkl
};

const char* label(const gemmi::Mtz::Column* col) {
  return col ? col->label.c_str() : "** N/A **";
}

DataForCheck read_refl_data(const char* input_path, const char* mapping, bool verbose) {
  using gemmi::Mtz;

  // parse mapping
  std::string mapping_str;
  const char* requested_cols[9] = {}; // fwt, phwt, delfwt, phdelwt, dfc, phdfc, fosc, fom, free
  if (mapping) {
    mapping_str = mapping;
    for (size_t start = 0, end = 0; end != std::string::npos; start = end + 1) {
      end = mapping_str.find(',', start);
      size_t colon = mapping_str.find(':', start);
      if (colon >= end)
        gemmi::fail("wrong argument for --check, missing ':'");
      mapping_str[colon] = '\0';
      if (end != std::string::npos)
        mapping_str[end] = '\0';
      size_t req_idx = 0;
      int key_id = 0;
      if (colon - start >= 3)
        key_id = gemmi::ialpha4_id(&mapping_str[start]);
      switch (key_id) {
        case gemmi::ialpha4_id("fwt"):  req_idx = 0; break;
        case gemmi::ialpha4_id("phwt"): req_idx = 1; break;
        case gemmi::ialpha4_id("delf"): req_idx = 2; break;
        case gemmi::ialpha4_id("phde"): req_idx = 3; break;
        case gemmi::ialpha4_id("dfc"):  req_idx = 4; break;
        case gemmi::ialpha4_id("phdf"): req_idx = 5; break;
        case gemmi::ialpha4_id("fosc"): req_idx = 6; break;
        case gemmi::ialpha4_id("fom"):  req_idx = 7; break;
        case gemmi::ialpha4_id("free"): req_idx = 8; break;
        default: gemmi::fail("Arg for --check has unknown key: ", start);
      }
      requested_cols[req_idx] = &mapping_str[colon + 1];
    }
  }

  DataForCheck result;
  if (verbose)
    printf("Reading reflections from %s ...\n", input_path);
  Mtz mtz;
  mtz.read_file_gz(input_path);
  gemmi::MtzDataProxy data_proxy{mtz};
  result.sg = data_proxy.spacegroup();
  if (!result.sg)
    gemmi::fail("unknown spacegroup in MTZ file");

  auto get_requested = [&](int idx) -> const Mtz::Column* {
    if (const char* name = requested_cols[idx])
      return &mtz.get_column_with_label(name);
    return nullptr;
  };

  printf("Columns used in checking map coefficients:\n");
  using FCols = std::array<const Mtz::Column*, 2>;
  auto get_columns = [&](const char* desc, int idx, std::initializer_list<const char*> labels) {
    FCols fcols = {get_requested(idx), get_requested(idx+1)};
    if (!fcols[0])
      fcols[0] = mtz.column_with_one_of_labels(labels, 'F');
    if (fcols[0] && !fcols[1])
      // Typically, the phase column follows directly the amplitude.
      // This assumption should be as robust as relying on column names.
      fcols[1] = fcols[0]->get_next_column_if_type('P');
    printf("    for %-20s %-10s %s\n", desc, label(fcols[0]), label(fcols[1]));
    return fcols;
  };
  // Using notation from Ian Tickle's mtzfix:
  // FM = 2-1 map coef,  FD = 1-1 map coef, FC = DFc
  FCols FM = get_columns("FM (normal map):", 0, {"FWT", "2FOFCWT"});
  FCols FD = get_columns("FD (difference map):", 2, {"DELFWT", "FOFCWT"});
  FCols FC = get_columns("D.Fc:", 4, {"DFC", "FC_ALL", "FC"});

  const Mtz::Column* Fo = get_requested(6);
  if (!Fo) {
    // Fobs is typically followed by sigma, not phase
    Fo = mtz.column_with_one_of_labels({"FOSC", "FP"}, 'F');
    if (Fo && Fo->get_next_column_if_type('P'))
      Fo = nullptr;
  }
  if (!Fo) {
    for (size_t i = 0; i < mtz.columns.size() - 1; ++i)
      if (mtz.columns[i].type == 'F' && mtz.columns[i+1].type == 'Q') {
        if (!Fo) {
          Fo = &mtz.columns[i];
        } else {
          Fo = nullptr;
          break;
        }
    }
  }
  printf("    for scaled Fo:                %s\n", label(Fo));
  result.has_fo = Fo != nullptr;

  const Mtz::Column* fom = get_requested(7);
  if (!fom)
    fom = mtz.column_with_one_of_labels({"FOM"}, 'W');
  printf("    for figure-of-merit m:        %s\n", label(fom));

  const Mtz::Column* free_flags = get_requested(8);
  if (!free_flags)
    free_flags = mtz.rfree_column();
  printf("    for free flags:               %s\n", label(free_flags));
  result.has_free_flags = free_flags != nullptr;

  if (!FM[1] || !FD[1] || !FC[1]) {
    printf("Checking aborted.\n");
    result.sg = nullptr;
    return result;
  }

  auto colidx = [](const Mtz::Column* col) { return col ? col->idx : 0; };
  std::array<size_t, 9> col_indices = {
    colidx(free_flags), // 0
    colidx(FM[0]),      // 1
    colidx(FM[1]),      // 2
    colidx(FD[0]),      // 3
    colidx(FD[1]),      // 4
    colidx(FC[0]),      // 5
    colidx(FC[1]),      // 6
    colidx(Fo),         // 7
    colidx(fom)         // 8
  };
  for (size_t i = 0; i < 9; ++i)
    if (col_indices[i] != 0)
      result.labels[i] = mtz.columns[col_indices[i]].label;

  result.refl.reserve(mtz.nreflections);

  auto colval = [&](size_t offset, size_t n, float missing=NAN) {
    return col_indices[n] != 0 ? data_proxy.get_num(offset + col_indices[n]) : missing;
  };
  for (size_t offset = 0; offset < data_proxy.size(); offset += data_proxy.stride()) {
    result.refl.push_back({
        data_proxy.get_hkl(offset),
        (int) colval(offset, 0, -1.f),
        colval(offset, 1),
        colval(offset, 2),
        colval(offset, 3),
        colval(offset, 4),
        colval(offset, 5),
        colval(offset, 6),
        colval(offset, 7),
        colval(offset, 8)});
  }
  return result;
}

void check_map_coef(DataForCheck& data) {
  gemmi::GroupOps gops = data.sg->operations();
  std::vector<int> zeroed_flags;

  if (data.has_free_flags) {
    std::unordered_map<int, unsigned> non_zero_count;
    for (ReflData& r : data.refl) {
      auto it = non_zero_count.emplace(r.free_flag, 0).first;
      if (r.fd != 0)  // neither 0 nor NaN
        ++it->second;
    }
    printf("Is FD (%s) set to 0 or NaN for any of %zu %s flags ...",
            data.labels[3].c_str(), non_zero_count.size(), data.labels[0].c_str());
    for (auto it : non_zero_count)
      if (it.second == 0)
        zeroed_flags.push_back(it.first);
    if (zeroed_flags.empty()) {
      printf(" no\n");
    } else {
      printf(" yes, for");
      for (int n : zeroed_flags)
        printf(" %d", n);
      printf("\n -> free reflections are NOT used for maps\n");
    }
  }

  // Check if phases of map coefficient and DFc are matching each other.
  // For opposite phases, negate the amplitude to make it comparable.
  printf("Phases %s and %s ... ", data.labels[2].c_str(), data.labels[4].c_str());
  for (ReflData& r : data.refl) {
    double phi_diff =  gemmi::angle_abs_diff(r.fm_phi, r.fd_phi, 360.);
    if (phi_diff > 0.5) {
      if (phi_diff > 180. - 0.5) {
        r.fd = -r.fd;  // make r.fd use the same phase as r.fm
      } else if (std::fabs(r.fm) > 1e-2 && std::fabs(r.fd) > 1e-2) {
        printf("DIFFER\nFor (%d %d %d) it is %g vs %g\n",
               r.hkl[0], r.hkl[1], r.hkl[2], r.fm_phi, r.fd_phi);
        printf("Why they differ? Aborting the check.\n");
        return;
      }
    }
  }
  printf("  match (mod 180)\n");
  printf("Phases %s and %s ... ", data.labels[6].c_str(), data.labels[2].c_str());
  for (ReflData& r : data.refl) {
    double phi_diff =  gemmi::angle_abs_diff(r.phi, r.fm_phi, 360.);
    if (phi_diff > 0.5) {
      if (phi_diff > 180. - 0.5) {
        r.fm = -r.fm;
        r.fd = -r.fd;
      } else if (std::fabs(r.fm) > 1e-2 && std::fabs(r.fc) > 1e-2) {
        printf("DIFFER\nFor (%d %d %d) it is %g vs %g\n",
                r.hkl[0], r.hkl[1], r.hkl[2], r.phi, r.fm_phi);
        printf("Probably Fc is not of a full model. Aborting the check.\n");
        return;
      }
    }
  }
  printf("  match (mod 180)\n");

  struct Counters {
    size_t count = 0;
    size_t fm1 = 0;  // FM != 2m.Fo - FC
    size_t fm2 = 0;  // FM != m.Fo
    size_t fm3 = 0;  // FM != FC
    size_t fmn = 0;  // FM is NaN
    size_t fd1 = 0;  // FD != FM - m.Fo
    size_t fd2 = 0;  // FD != 2(FM - m.Fo)
    size_t fd3 = 0;  // FD != m.Fo - FC
    size_t fd4 = 0;  // FD != 2(m.Fo - FC)
    size_t fd5 = 0;  // FD != FM - FC
    size_t fd6 = 0;  // FD != (FM - FC)/2
    size_t fd0 = 0;  // FD != 0
    size_t fdn = 0;  // FD is NaN

    std::string get_fm() {
      if (fmn == count)
        return " = NaN";
      std::string s;
      if (fm2 == 0) s += " = m.Fo";
      if (fm3 == 0) s += " = D.Fc";
      // if FM = m.Fo = D.Fc there is no need to add "= 2m.Fo - D.Fc"
      else if (fm1 == 0) s += " = 2m.Fo - D.Fc";
      if (s.empty())
          s = " = ?";
      if (fmn != 0)
        s += " (if not NaN)";
      return s;
    }

    std::string get_fd() {
      if (fdn == count)
        return " = NaN";
      std::string s;
      if (fd0 == 0) {
        s = " = 0";
      } else {
        if (fd3 == 0) s += " = m.Fo - D.Fc";
        if (fd4 == 0) s += " = 2(m.Fo - D.Fc)";
        if (fd1 == 0) s += " = FM - m.Fo";
        if (fd2 == 0) s += " = 2(FM - m.Fo)";
        if (fd5 == 0) s += " = FM - D.Fc";
        if (fd6 == 0) s += " = (FM - D.Fc)/2";
      }
      if (s.empty())
        s = " = ?";
      if (fdn != 0)
        s += " (if not NaN)";
      return s;
    }

    void add_refl(const ReflData& r) {
      ++count;
      float mFo = r.m * r.fo;
      if (std::isnan(r.fm)) {
        fmn++;
      } else {
        fm1 += !(std::fabs(r.fm - (2 * mFo - r.fc)) < 0.1f);
        fm2 += !(std::fabs(r.fm - mFo)              < 0.1f);
        fm3 += !(std::fabs(r.fm - r.fc)             < 0.01f);
      }
      if (std::isnan(r.fd)) {
        fdn++;
      } else {
        fd1 += !(std::fabs(r.fd - (r.fm - mFo))     < 0.1f);
        fd2 += !(std::fabs(r.fd - 2 * (r.fm - mFo)) < 0.1f);
        fd3 += !(std::fabs(r.fd - (mFo - r.fc))     < 0.1f);
        fd4 += !(std::fabs(r.fd - 2 * (mFo - r.fc)) < 0.1f);
        fd5 += !(std::fabs(r.fd - (r.fm - r.fc))    < 0.1f);
        fd6 += !(std::fabs(r.fd - 0.5f*(r.fm - r.fc)) < 0.1f);
        fd0 += !(std::fabs(r.fd - 0)                < 0.01f);
      }
    }
  };

  Counters counters[3]; // unused, centric, acentric

  for (ReflData& r : data.refl) {
    int idx = 2;
    if (r.m == 0 || (data.has_fo && std::isnan(r.fo)) ||
        gemmi::in_vector(r.free_flag, zeroed_flags))
      idx = 0;
    else if (gops.is_reflection_centric(r.hkl))
      idx = 1;
    counters[idx].add_refl(r);
  }
  const char* desc[3] = {
    "missing/unused reflections",
    "centric reflections (excl. missing/unused)",
    "acentric reflections (excl. missing/unused)"
  };

  for (int i = 2; i >= 0; --i) {
    Counters& c = counters[i];
    printf("For all %zu %s:\n", c.count, desc[i]);
    if (c.count == 0) {
      printf("    n/a\n");
    } else {
      printf("    FM%s\n", c.get_fm().c_str());
      printf("    FD%s\n", c.get_fd().c_str());
    }
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.options[Check]) {
    if (p.nonOptionsCount() != 1) {
      fprintf(stderr, "%s --check requires 1 arguments, got %d.",
              p.program_name, p.nonOptionsCount());
      p.print_try_help_and_exit("");
    }
  } else if (p.options[GridQuery]) {  // 1 arg used, or 2 args if combined with other options
    p.require_input_files_as_args();
  } else {
    p.require_positional_args(2);
  }
  try {
    if (p.options[Check]) {
      bool verbose = p.options[Verbose];
      DataForCheck data = read_refl_data(p.nonOption(0), p.options[Check].arg, verbose);
      if (data.sg)
        check_map_coef(data);
      // --check and -G could be combined
      if (!p.options[GridQuery])
        return 0;
    }
    transform_sf_to_map(p);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
