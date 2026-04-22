// Copyright Global Phasing Ltd.

/// @file
/// @brief DDL1/DDL2 dictionary-based validation of CIF and mmCIF files.

#ifndef GEMMI_DDL_HPP_
#define GEMMI_DDL_HPP_

#include <map>
#include <memory>  // for unique_ptr
#include <regex>
#include "cifdoc.hpp"  // for cif::Document
#include "logger.hpp"  // for Logger

namespace gemmi { namespace cif {

/// Represents a CIF dictionary (DDL1 or DDL2 ontology) for validation.
///
/// A DDL (Data Definition Language) dictionary defines the structure, constraints,
/// and validation rules for CIF data. This class can load and use either:
/// - **DDL1** dictionaries (IUCr core, chemical structures)
/// - **DDL2** dictionaries (macromolecular CIF / mmCIF, used by PDB)
///
/// After loading a dictionary with read_ddl(), you can validate CIF documents
/// against it to check for missing mandatory items, type violations, enumeration
/// violations, unique key violations, and other data integrity issues.
struct GEMMI_DLL Ddl {
  /// Logger for validation messages and warnings.
  /// Member functions use this logger's callback and threshold settings for output.
  Logger logger;

  // Configuration flags - set these before calling read_ddl()

  /// Report unknown tags (tags not defined in the dictionary).
  /// Useful for catching typos in tag names.
  bool print_unknown_tags = true;

  // The following flags apply to DDL2 dictionaries only

  /// Enable validation using regular expression patterns (DDL2 _item_type.code).
  bool use_regex = true;

  /// Use context-dependent validation rules (DDL2).
  /// If true, validates items in specific category contexts.
  bool use_context = false;

  /// Use parent-child item relationships (DDL2 _item_linked).
  /// If true, enforces dependencies between items.
  bool use_parents = false;

  /// Validate mandatory items (DDL2 _item.mandatory_code).
  /// If true, reports missing items marked as mandatory.
  bool use_mandatory = true;

  /// Validate unique keys in loops (DDL2 _item_linked.key_id).
  /// If true, checks that unique key values don't repeat.
  bool use_unique_keys = true;

  /// Use PDBx deposition-specific validation checks.
  /// If true, uses _pdbx-prefixed dictionary items instead of standard ones
  /// (_pdbx_item_type.code instead of _item_type.code, etc.).
  /// This mode is typically used during structure deposition to PDB.
  bool use_deposition_checks = false;

  // Read-only fields set when reading a dictionary

  /// Major version of the loaded DDL (1 or 2).
  /// Read from _dictionary_version or similar field in the dictionary.
  int major_version = 0;

  /// Name of the dictionary (e.g., "cif_core.dic" or "mmcif_pdbx_v50").
  /// Read from _dictionary_name (DDL1) or _dictionary.title (DDL2).
  std::string dict_name;

  /// Version string of the dictionary (e.g., "2.0.11").
  /// Read from _dictionary_version or _dictionary.version.
  std::string dict_version;

  Ddl() = default;

  // Copy/assignment deleted: MSVC dllexport cannot handle the unique_ptr
  // member in ddl_docs_. Instances should be moved or held in stable storage.
  Ddl(Ddl const&) = delete;
  Ddl& operator=(Ddl const&) = delete;

  /// Load a DDL dictionary into this validator.
  ///
  /// Parses a DDL1 or DDL2 dictionary document and indexes it for validation.
  /// The document is moved into internal storage to manage its lifetime.
  /// The dictionary version (DDL1 or DDL2) is auto-detected.
  ///
  /// Configuration flags (e.g., use_mandatory, use_regex) should be set
  /// before calling this function.
  ///
  /// @param doc CIF document containing a DDL dictionary (will be moved)
  void read_ddl(cif::Document&& doc);

  /// Validate all blocks in a CIF document against this dictionary.
  ///
  /// Checks all blocks in the document and reports validation errors
  /// via the configured logger.
  ///
  /// @param doc The CIF document to validate
  /// @return true if validation passes, false if errors are found
  ///
  /// @see validate_block() to validate individual blocks
  bool validate_cif(const cif::Document& doc) const;

  /// Validate a single CIF block against this dictionary.
  ///
  /// Performs all enabled validation checks on the block:
  /// - Mandatory items (if use_mandatory=true)
  /// - Item types and enumeration values
  /// - Regular expression patterns (if use_regex=true)
  /// - Unique keys (if use_unique_keys=true)
  /// - Parent-child relationships (if use_parents=true)
  /// - Unknown tags (if print_unknown_tags=true)
  ///
  /// @param b      The CIF block to validate
  /// @param source Source identifier for error messages (e.g., block name or filename)
  /// @return true if validation passes, false if errors are found
  bool validate_block(const cif::Block& b, const std::string& source) const;

  /// Check audit conformance fields in a CIF document.
  ///
  /// Verifies that the document's audit records match dictionary expectations
  /// (e.g., _audit_conform_dict_name, _audit_conform_dict_version).
  /// Reports mismatches via the logger.
  ///
  /// @param doc The CIF document to check
  void check_audit_conform(const cif::Document& doc) const;

  /// Access the regex patterns loaded from the dictionary.
  ///
  /// Returns a map of tag names to compiled regular expressions
  /// that constrain the format of values for those tags (DDL2 validation).
  ///
  /// @return Map of regex patterns indexed by tag name
  const std::map<std::string, std::regex>& regexes() const { return regexes_; }

private:
  /// Internal representation of DDL2 parent-child item relationships.
  ///
  /// Links parent and child tags that must be coordinated in the data.
  /// Used for enforcing referential integrity (use_parents=true).
  struct ParentLink {
    std::string group;                   ///< Name of the linked group
    std::vector<std::string> child_tags; ///< Child item tags
    std::vector<std::string> parent_tags;///< Parent item tags
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
