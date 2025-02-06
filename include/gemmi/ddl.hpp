// Copyright Global Phasing Ltd.
//
// Using DDL1/DDL2 dictionaries to validate CIF/mmCIF files.

#ifndef GEMMI_DDL_HPP_
#define GEMMI_DDL_HPP_

#include <map>
#include <memory>  // for unique_ptr
#include <regex>
#include "cifdoc.hpp"  // for cif::Document
#include "logger.hpp"  // for Logger

namespace gemmi { namespace cif {

/// Represents DDL1 or DDL2 dictionary (ontology).
struct GEMMI_DLL Ddl {
  /// member functions use logger's callback and threshold for output
  Logger logger;
  // configuration - some of these flag must be set before read_ddl()
  bool print_unknown_tags = true;
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

  /// it moves doc to ddl_docs_ to control lifetime and prevent modifications
  void read_ddl(cif::Document&& doc);

  bool validate_cif(const cif::Document& doc) const;
  bool validate_block(const cif::Block& b, const std::string& source) const;

  void check_audit_conform(const cif::Document& doc) const;

  const std::map<std::string, std::regex>& regexes() const { return regexes_; }

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
  void check_mandatory_items(const cif::Block& b) const;
  void check_unique_keys_in_loop(const cif::Loop& loop, const cif::Block& block) const;
  void check_parents(const cif::Block& b) const;
  void check_parent_link(const ParentLink& link, const cif::Block& b) const;
  void read_ddl1_block(cif::Block& block);
  void read_ddl2_block(cif::Block& block);

  template<class... Args> void warn(const cif::Block& b, Args const&... args) const {
    logger.level<3>('[', b.name, "] ", args...);
  }
};

}} // namespace gemmi::cif
#endif
