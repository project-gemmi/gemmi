// Copyright Global Phasing Ltd.
//
// Using DDL1/DDL2 dictionaries to validate CIF/mmCIF files.

#ifndef GEMMI_DDL_HPP_
#define GEMMI_DDL_HPP_

#include <map>
#include <memory>  // for unique_ptr
#include <ostream>
#include <regex>
#include "cifdoc.hpp"  // for cif::Document

namespace gemmi { namespace cif {

/// Represents DDL1 or DDL2 dictionary (ontology).
struct GEMMI_DLL Ddl {
  // configuration - some of these flag must be set before read_ddl()
  bool print_unknown_tags = true;
  bool print_extra_diagnostics = false;
  // these flags below are relevant to DDL2 only
  bool use_regex = true;
  bool use_context = false;
  bool use_parents = false;
  bool use_mandatory = true;
  bool use_unique_keys = true;

  // variables set when reading DLL; normally, no need to change them
  int major_version = 0;  // currently 1 and 2 are supported
  std::string dict_name;  // _dictionary_name or _dictionary.title
  std::string dict_version;  // _dictionary_version or _dictionary.version

  Ddl() = default;
  // MSVC with dllexport attempts to export all non-deleted member functions,
  // failing with Error C2280 (because of ddl_docs_) if we don't delete these:
  Ddl(Ddl const&) = delete;
  Ddl& operator=(Ddl const&) = delete;

  void read_ddl(cif::Document&& doc, std::ostream& out);

  bool validate_cif(const cif::Document& doc, std::ostream& out) const;

  void check_audit_conform(const cif::Document& doc, std::ostream& out) const;

private:
  // items from DDL2 _pdbx_item_linked_group[_list]
  struct ParentLink {
    std::string group;
    std::vector<std::string> child_tags;
    std::vector<std::string> parent_tags;
  };

  std::vector<std::unique_ptr<cif::Document>> ddl_docs_;
  std::map<std::string, cif::Block*> name_index_;
  std::map<std::string, std::regex> regexes_;
  std::vector<ParentLink> parents_;
  // storage for DDL2 _item_linked.child_name -> _item_linked.parent_name
  std::map<std::string, std::string> item_parents_;

  cif::Block* find_rules(const std::string& name) const {
    auto iter = name_index_.find(to_lower(name));
    return iter != name_index_.end() ? iter->second : nullptr;
  }
  void check_mandatory_items(const cif::Block& b, std::ostream& out) const;
  void check_unique_keys_in_loop(const cif::Loop& loop, std::ostream& out,
                                 const std::string& block_name) const;
  void check_parents(const cif::Block& b, std::ostream& out) const;
  void read_ddl1_block(cif::Block& block);
  void read_ddl2_block(cif::Block& block, std::ostream& out);
};

}} // namespace gemmi::cif
#endif
