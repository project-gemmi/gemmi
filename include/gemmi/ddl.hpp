//! @file
//! @brief Using DDL1/DDL2 dictionaries to validate CIF/mmCIF files.

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

//! @brief Represents DDL1 or DDL2 dictionary (ontology).
//!
//! DDL (Dictionary Definition Language) dictionaries define the structure
//! and constraints for CIF/mmCIF files. This class loads and applies these
//! validation rules.
struct GEMMI_DLL Ddl {
  Logger logger;  //!< Member functions use logger's callback and threshold for output

  // Configuration - some of these flags must be set before read_ddl()
  bool print_unknown_tags = true;  //!< Print warnings for unknown tags

  // Flags relevant to DDL2 only
  bool use_regex = true;           //!< Use regex patterns for validation
  bool use_context = false;        //!< Use context-dependent validation
  bool use_parents = false;        //!< Check parent-child relationships
  bool use_mandatory = true;       //!< Check mandatory items
  bool use_unique_keys = true;     //!< Check unique key constraints
  bool use_deposition_checks = false;  //!< Use _pdbx-prefixed validation items

  // Variables set when reading DDL; normally, no need to change them
  int major_version = 0;        //!< DDL version (1 or 2)
  std::string dict_name;        //!< _dictionary_name or _dictionary.title
  std::string dict_version;     //!< _dictionary_version or _dictionary.version

  Ddl() = default;
  // MSVC with dllexport attempts to export all non-deleted member functions,
  // failing with Error C2280 (because of ddl_docs_) if we don't delete these:
  Ddl(Ddl const&) = delete;
  Ddl& operator=(Ddl const&) = delete;

  //! @brief Read DDL dictionary.
  //! @param doc DDL dictionary document (moved to control lifetime)
  //!
  //! It moves doc to ddl_docs_ to control lifetime and prevent modifications.
  void read_ddl(cif::Document&& doc);

  //! @brief Validate CIF document against dictionary.
  //! @param doc CIF document to validate
  //! @return True if valid
  bool validate_cif(const cif::Document& doc) const;

  //! @brief Validate single CIF block.
  //! @param b Block to validate
  //! @param source Source description for error messages
  //! @return True if valid
  bool validate_block(const cif::Block& b, const std::string& source) const;

  //! @brief Check _audit_conform entries.
  //! @param doc CIF document to check
  void check_audit_conform(const cif::Document& doc) const;

  //! @brief Get compiled regex patterns.
  //! @return Map of tag names to regex patterns
  const std::map<std::string, std::regex>& regexes() const { return regexes_; }

private:
  //! @brief Items from DDL2 _pdbx_item_linked_group[_list].
  //!
  //! Represents parent-child relationships between tags.
  struct ParentLink {
    std::string group;                        //!< Link group name
    std::vector<std::string> child_tags;      //!< Child tag names
    std::vector<std::string> parent_tags;     //!< Parent tag names
  };

  std::vector<std::unique_ptr<cif::Document>> ddl_docs_;  //!< Owned DDL documents
  std::map<std::string, cif::Block*> name_index_;         //!< Tag name to block lookup
  std::map<std::string, std::regex> regexes_;             //!< Compiled regex patterns
  std::vector<ParentLink> parents_;                       //!< Parent-child links
  std::map<std::string, std::string> item_parents_;       //!< DDL2 _item_linked mappings

  //! @brief Find validation rules for tag.
  //! @param name Tag name
  //! @return Block containing rules or nullptr
  cif::Block* find_rules(const std::string& name) const {
    auto iter = name_index_.find(to_lower(name));
    return iter != name_index_.end() ? iter->second : nullptr;
  }

  //! @brief Check for mandatory items in block.
  //! @param b Block to check
  void check_mandatory_items(const cif::Block& b) const;

  //! @brief Check unique key constraints in loop.
  //! @param loop Loop to check
  //! @param block Containing block
  void check_unique_keys_in_loop(const cif::Loop& loop, const cif::Block& block) const;

  //! @brief Check parent-child relationships.
  //! @param b Block to check
  void check_parents(const cif::Block& b) const;

  //! @brief Check specific parent link.
  //! @param link Parent link to check
  //! @param b Block to check
  void check_parent_link(const ParentLink& link, const cif::Block& b) const;

  //! @brief Read DDL1 dictionary block.
  //! @param block Block to read
  void read_ddl1_block(cif::Block& block);

  //! @brief Read DDL2 dictionary block.
  //! @param block Block to read
  void read_ddl2_block(cif::Block& block);

  //! @brief Log warning message.
  //! @tparam Args Argument types
  //! @param b Block being validated
  //! @param args Message parts
  template<class... Args> void warn(const cif::Block& b, Args const&... args) const {
    logger.level<3>('[', b.name, "] ", args...);
  }
};

}} // namespace gemmi::cif
#endif
