// Copyright 2017-2020 Global Phasing Ltd.

#include "gemmi/ddl.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/numb.hpp"
#include "gemmi/dirwalk.hpp"  // for CifWalk
#include "gemmi/read_cif.hpp"  // for read_cif_gz, check_cif_syntax_gz
#include "validate_mon.h"  // for check_monomer_doc
#include <stdio.h>
#include <stdexcept>  // for std::runtime_error

#ifdef GEMMI_ANALYZE_RULES
# include "gemmi/cif.hpp"
# include <tao/pegtl/analyze.hpp>
#endif

#define GEMMI_PROG validate
#include "options.h"

namespace cif = gemmi::cif;

namespace {

enum OptionIndex {
  Quiet=4, Fast, Stat, Context, Ddl, NoRegex, NoMandatory, NoUniqueKeys,
  Parents, Recurse, Monomer, Zscore, Ccd, AuditDate
};
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None, "Usage: " EXE_NAME " [options] FILE [...]"
                                "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Quiet, 0, "q", "quiet", Arg::None, "  -q, --quiet  \tShow only errors." },
  { Fast, 0, "f", "fast", Arg::None, "  -f, --fast  \tSyntax-only check." },
  { Stat, 0, "s", "stat", Arg::None, "  -s, --stat  \tShow token statistics" },
  { Recurse, 0, "r", "recursive", Arg::None,
    "  -r, --recursive  \tRecurse directories and process all CIF files." },
  { Ddl, 0, "d", "ddl", Arg::Required, "  -d, --ddl=PATH  \tDDL for validation." },

  { NoOp, 0, "", "", Arg::None, "\nOptional checks:" },
  { Context, 0, "c", "context", Arg::None,
    "  -c, --context  \tCheck _pdbx_{category|item}_context.type." },
  { NoRegex, 0, "", "no-regex", Arg::None,
    "  --no-regex  \tSkip regex checking (when using DDL2)" },
  { NoMandatory, 0, "", "no-mandatory", Arg::None,
    "  --no-mandatory  \tSkip checking if mandatory tags are present." },
  { NoUniqueKeys, 0, "", "no-unique", Arg::None,
    "  --no-unique  \tSkip checking if category keys (DDL2) are unique." },
  { Parents, 0, "p", "", Arg::None,
    "  -p  \tCheck if parent items (DDL2) are present." },

  { NoOp, 0, "", "", Arg::None, "\nValidation specific to CCP4 monomer files:" },
  { Monomer, 0, "m", "monomer", Arg::None,
    "  -m, --monomer  \tRun checks specific to monomer dictionary." },
  { Zscore, 0, "", "z-score", Arg::Float,
    "  --z-score=Z  \tUse Z for validating _chem_comp_atom.[xyz] (default: 2.0)." },
  { Ccd, 0, "", "ccd", Arg::Required,
    "  --ccd=PATH  \tCCD file for comparison." },
  { AuditDate, 0, "", "audit-on", Arg::Required,
    "  --audit-on=DATE  \tCheck only if CCD component was remediated on DATE." },
  { 0, 0, 0, 0, 0, 0 }
};

// basic types, used for token statistics only
enum class ValueType : unsigned char {
  NotSet,
  Char,
  Numb,
  Dot,
  QuestionMark,
};

inline std::string value_type_to_str(ValueType v) {
  switch (v) {
    case ValueType::NotSet: return "n/a";
    case ValueType::Char: return "char";
    case ValueType::Numb: return "numb";
    case ValueType::Dot: return "'.'";
    case ValueType::QuestionMark: return "'?'";
  }
  return "";
}

// For now the infer_* functions are used only here, not sure where they belong
inline ValueType infer_value_type(const std::string& val) {
  assert(!val.empty());
  if (val == ".")
    return ValueType::Dot;
  if (val == "?")
    return ValueType::QuestionMark;
  if (cif::is_numb(val))
    return ValueType::Numb;
  return ValueType::Char;
}

std::string format_7zd(size_t k) {
  char buf[64];
  snprintf(buf, 63, "%7zu", k);
  return buf;
}

std::string token_stats(const cif::Document& d) {
  size_t nframes = 0, nvals = 0, nloops = 0, nlooptags = 0, nloopvals = 0;
  size_t vals_by_type[5] = {0};
  size_t looptags_by_type[5] = {0};
  for (const cif::Block& block : d.blocks) {
    for (const cif::Item& item : block.items) {
      if (item.type == cif::ItemType::Pair) {
        nvals++;
        ValueType vt = infer_value_type(item.pair[1]);
        vals_by_type[static_cast<int>(vt)]++;
      } else if (item.type == cif::ItemType::Frame) {
        nframes++;
      } else if (item.type == cif::ItemType::Loop) {
        nloops++;
        size_t width = item.loop.width();
        nlooptags += width;
        nloopvals += item.loop.values.size();
        for (size_t i = 0; i != width; ++i) {
          ValueType vt = ValueType::NotSet;
          // TODO: ConstColumn(const::Item*, ...)
          const cif::Column col(const_cast<cif::Item*>(&item), i);
          for (const std::string& v : col) {
            ValueType this_vt = infer_value_type(v);
            if (this_vt != vt) {
              // if we are here: vt != ValueType::Char
              if (vt == ValueType::NotSet || this_vt == ValueType::Numb) {
                vt = this_vt;
              } else if (this_vt == ValueType::Char) {
                vt = this_vt;
                break;
              }
            }
          }
          looptags_by_type[static_cast<int>(vt)]++;
        }
      }
    }
  }
  std::string info;
  gemmi::cat_to(info, format_7zd(d.blocks.size()), " block(s)\n");
  gemmi::cat_to(info, format_7zd(nframes), " frames\n");
  gemmi::cat_to(info, format_7zd(nvals), " non-loop items:");
  for (int i = 1; i != 5; ++i)
    gemmi::cat_to(info, "  ", value_type_to_str(static_cast<ValueType>(i)),
                  ':', vals_by_type[i]);
  gemmi::cat_to(info, '\n', format_7zd(nloops), " loops w/"
                "\n        ", format_7zd(nlooptags), " tags:");
  for (int i = 1; i != 5; ++i)
    gemmi::cat_to(info, "  ", value_type_to_str(static_cast<ValueType>(i)),
                  ':', looptags_by_type[i]);
  gemmi::cat_to(info, "\n        ", format_7zd(nloopvals), " values\n");
  return info;
}

// Empty loop is not a valid CIF syntax, but we parse it to accommodate
// some broken CIF files. Only validation shows an error.
void check_empty_loops(const cif::Block& block) {
  for (const cif::Item& item : block.items) {
    if (item.type == cif::ItemType::Loop) {
      if (item.loop.values.empty() && !item.loop.tags.empty())
        throw std::runtime_error("Empty loop in block " + block.name +
                                 ": " + item.loop.tags[0]);
    } else if (item.type == cif::ItemType::Frame) {
      check_empty_loops(item.frame);
    }
  }
}

bool column_contains(const cif::Column& column, const char* value) {
  for (const std::string& v : column)
    if (cif::as_string(v) == value)
      return true;
  return false;
}

bool process_file(const char* path, const cif::Ddl& dict,
                  const std::map<std::string, cif::Block>& ccd_map,
                  const std::vector<option::Option>& options) {
  bool ok = true;
  std::string msg;
  if (options[Verbose])
    printf("Reading %s...\n", path);
  try {
    if (options[Fast]) {
      ok = gemmi::check_cif_syntax_gz(path, &msg);
    } else {
      cif::Document doc = gemmi::read_cif_gz(path);
      for (const cif::Block& block : doc.blocks) {
        if (block.name == " ")
          printf("%s: missing block name (bare data_)\n", doc.source.c_str());
        check_empty_loops(block);
      }
      if (options[Stat])
        msg = token_stats(doc);
      if (options[Ddl]) {
        dict.check_audit_conform(doc);
        ok = dict.validate_cif(doc);
      }

      if (options[Monomer] || options[Zscore] || options[Ccd]) {
        for (const cif::Block& block : doc.blocks) {
          if (block.name == "comp_list")
            continue;
          const std::string cc_name =
            block.name.substr(gemmi::starts_with(block.name, "comp_") ? 5 : 0);
          try {
            if (options[Monomer] || options[Zscore]) {
              double z_score = options[Zscore] ? std::atof(options[Zscore].arg) : 2.0;
              check_monomer(block, z_score);
            }
            if (options[Ccd]) {
              auto it = ccd_map.find(cc_name);
              if (it == ccd_map.end()) {
                printf("%s [ccd] monomer not found in CCD file(s)\n", cc_name.c_str());
                continue;
              }
              const cif::Block& ccd_block = it->second;
              if (options[AuditDate]) {
                cif::Column col = const_cast<cif::Block&>(ccd_block)
                                   .find_values("_pdbx_chem_comp_audit.date");
                if (!column_contains(col, options[AuditDate].arg)) {
                  if (options[Verbose])
                    printf("%s [ccd] ignored - audit date does not match\n", cc_name.c_str());
                  continue;
                }
              }
              compare_monomer_with_ccd(block, ccd_block, options[Verbose]);
            }
          } catch (const std::exception& e) {
            fprintf(stderr, "Failed to interpret monomer: block %s from %s\n%s\n",
                    block.name.c_str(), doc.source.c_str(), e.what());
          }
        }
      }
    }
  } catch (std::runtime_error& e) {
    ok = false;
    msg = e.what();
  }
  if (!msg.empty())
    printf("%s\n", msg.c_str());

  if (options[Verbose])
    puts(ok ? "OK" : "FAILED");
  return ok;
}

} // anonymous namespace


int GEMMI_MAIN(int argc, char **argv) {
#ifdef GEMMI_ANALYZE_RULES // for debugging only
  tao::pegtl::analyze<cif::rules::file>();
#endif
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();

  bool total_ok = true;
  cif::Ddl dict;
  dict.logger = {&gemmi::Logger::to_stdout, 5 + p.options[Verbose].count()};
  dict.print_unknown_tags = !p.options[Quiet];
  dict.use_regex = !p.options[NoRegex];
  dict.use_context = p.options[Context];
  dict.use_parents = p.options[Parents];
  dict.use_mandatory = !p.options[NoMandatory];
  dict.use_unique_keys = !p.options[NoUniqueKeys];
  if (p.options[Ddl]) {
    try {
      for (option::Option* ddl = p.options[Ddl]; ddl; ddl = ddl->next())
        dict.read_ddl(gemmi::read_cif_gz(ddl->arg));
    } catch (std::runtime_error& e) {
      fprintf(stderr, "Error when reading dictionary: %s\n", e.what());
      return EXIT_FAILURE;
    }
  }
  std::map<std::string, cif::Block> ccd_map;
  if (p.options[Ccd]) {
    try {
      for (option::Option* ccd = p.options[Ccd]; ccd; ccd = ccd->next()) {
        cif::Document ccd_doc(gemmi::read_cif_gz(ccd->arg));
        for (cif::Block& ccd_block : ccd_doc.blocks) {
          std::string name = ccd_block.name;
          ccd_map.emplace(name, std::move(ccd_block));
        }
      }
    } catch (std::runtime_error& e) {
      fprintf(stderr, "Error when reading CCD file: %s\n", e.what());
      return EXIT_FAILURE;
    }
  }
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* path = p.nonOption(i);
    if (p.options[Recurse]) {
      for (const std::string& file : gemmi::CifWalk(path))
        total_ok = process_file(file.c_str(), dict, ccd_map, p.options) && total_ok;
    } else {
      total_ok = process_file(path, dict, ccd_map, p.options) && total_ok;
    }
  }
  return total_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
