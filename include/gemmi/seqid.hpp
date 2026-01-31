//! @file
//! @brief Sequence identifiers (residue numbers with insertion codes).
//!
//! SeqId -- residue number and insertion code together.

// Copyright 2017 Global Phasing Ltd.
//
// SeqId -- residue number and insertion code together.

#ifndef GEMMI_SEQID_HPP_
#define GEMMI_SEQID_HPP_

#include <climits>    // for INT_MIN
#include <cstdlib>    // for strtol
#include <stdexcept>  // for invalid_argument
#include <string>
#include "util.hpp"   // for cat

namespace gemmi {

//! @brief Optional integer with special "not-set" value.
//! @tparam N The value representing "not-set"
//!
//! Used for sequence numbers that may be absent.
template<int N> struct OptionalInt {
  enum { None=N };  //!< Special value indicating "not-set"
  int value = None;  //!< Stored value

  OptionalInt() = default;
  OptionalInt(int n) : value(n) {}
  bool has_value() const { return value != None; }
  std::string str(char null='?') const {
    return has_value() ? std::to_string(value) : std::string(1, null);
  }
  OptionalInt& operator=(int n) { value = n; return *this; }
  bool operator==(const OptionalInt& o) const { return value == o.value; }
  bool operator!=(const OptionalInt& o) const { return value != o.value; }
  bool operator<(const OptionalInt& o) const {
    return has_value() && o.has_value() && value < o.value;
  }
  bool operator==(int n) const { return value == n; }
  bool operator!=(int n) const { return value != n; }
  OptionalInt operator+(OptionalInt o) const {
    return OptionalInt(has_value() && o.has_value() ? value + o.value : N);
  }
  OptionalInt operator-(OptionalInt o) const {
    return OptionalInt(has_value() && o.has_value() ? value - o.value : N);
  }
  OptionalInt& operator+=(int n) { if (has_value()) value += n; return *this; }
  OptionalInt& operator-=(int n) { if (has_value()) value -= n; return *this; }
  explicit operator int() const { return value; }
  explicit operator bool() const { return has_value(); }
  // these are defined for partial compatibility with C++17 std::optional
  using value_type = int;
  int& operator*() { return value; }
  const int& operator*() const { return value; }
  int& emplace(int n) { value = n; return value; }
  void reset() noexcept { value = None; }
};

//! @brief Sequence identifier (residue number + insertion code).
//!
//! Uniquely identifies a residue position in a PDB file.
//! Follows PDB format conventions for residue numbering.
struct SeqId {
  using OptionalNum = OptionalInt<INT_MIN>;

  OptionalNum num;   //!< Sequence number
  char icode = ' ';  //!< Insertion code (space if none)

  SeqId() = default;
  SeqId(int num_, char icode_) { num = num_; icode = icode_; }
  SeqId(OptionalNum num_, char icode_) { num = num_; icode = icode_; }

  //! @brief Construct from string representation (e.g., "123", "45A").
  //! @param str String containing number and optional insertion code
  //! @throws std::invalid_argument if string is not valid
  explicit SeqId(const std::string& str) {
    char* endptr;
    num = std::strtol(str.c_str(), &endptr, 10);
    if (endptr == str.c_str() || (*endptr != '\0' && endptr[1] != '\0'))
      throw std::invalid_argument("Not a seqid: " + str);
    icode = (*endptr | 0x20);
  }

  bool operator==(const SeqId& o) const {
    return num == o.num && ((icode ^ o.icode) & ~0x20) == 0;
  }
  bool operator!=(const SeqId& o) const { return !operator==(o); }
  bool operator<(const SeqId& o) const {
    return (*num * 256 + icode) < (*o.num * 256 + o.icode);
  }
  bool operator<=(const SeqId& o) const { return !(o < *this); }

  char has_icode() const { return icode != ' '; }

  std::string str(bool dot_before_icode=false) const {
    std::string r = num.str();
    if (icode != ' ') {
      if (dot_before_icode)
        r += '.';
      r += icode;
    }
    return r;
  }
};

//! @brief Complete residue identifier.
//!
//! Combines sequence ID, residue name, and segment ID.
//! Uniquely identifies a residue in a chain.
struct ResidueId {
  SeqId seqid;  //!< Sequence number + insertion code
  std::string segment;  //!< Segment ID (up to 4 chars in PDB)
  std::string name;  //!< Residue name (e.g., "ALA", "GLY")

  //! @brief Get grouping key for residue.
  //! @return SeqId (used for first_conformation iterators, etc.)
  //!
  //! used for first_conformation iterators, etc.
  SeqId group_key() const { return seqid; }

  //! @brief Check if residue IDs match (including segment).
  //! @param o Other residue ID
  //! @return True if all fields match
  bool matches(const ResidueId& o) const {
    return seqid == o.seqid && segment == o.segment && name == o.name;
  }
  //! @brief Check if residue IDs match (ignoring segment).
  //! @param o Other residue ID
  //! @return True if seqid and name match
  bool matches_noseg(const ResidueId& o) const {
    return seqid == o.seqid && name == o.name;
  }
  bool operator==(const ResidueId& o) const { return matches(o); }
  //! @brief Convert to string representation.
  //! @return String like "123(ALA)"
  std::string str() const { return cat(seqid.str(), '(', name, ')'); }
};

//! @brief Format atom identifier string.
//! @param chain_name Chain identifier
//! @param res_id Residue identifier
//! @param atom_name Atom name
//! @param altloc Alternate location indicator
//! @param as_cid Format as CID (mmdb selection syntax)
//! @return Formatted atom string
inline std::string atom_str(const std::string& chain_name,
                            const ResidueId& res_id,
                            const std::string& atom_name,
                            char altloc,
                            bool as_cid=false) {
  std::string r = as_cid ? "//" + chain_name : chain_name;
  r += '/';
  if (!as_cid) {
    r += res_id.name;
    r += ' ';
  }
  r += res_id.seqid.str(as_cid);
  if (as_cid && atom_name == "null")
    return r;
  r += '/';
  r += atom_name;
  if (altloc) {
    r += as_cid ? ':' : '.';
    r += altloc;
  }
  return r;
}

//! @brief Complete atom address (chain + residue + atom + altloc).
//!
//! Uniquely identifies a specific atom in a structure.
//! Includes all hierarchical information needed to locate an atom.
struct AtomAddress {
  std::string chain_name;  //!< Chain identifier
  ResidueId res_id;  //!< Residue identifier
  std::string atom_name;  //!< Atom name (e.g., "CA", "N")
  char altloc = '\0';  //!< Alternate location indicator

  AtomAddress() = default;
  AtomAddress(const std::string& ch, const ResidueId& resid,
              const std::string& atom, char alt='\0')
    : chain_name(ch), res_id(resid), atom_name(atom), altloc(alt) {}
  AtomAddress(const std::string& ch, const SeqId& seqid, const std::string& res,
              const std::string& atom, char alt='\0')
    : chain_name(ch), res_id({seqid, "", res}), atom_name(atom), altloc(alt) {}
  bool operator==(const AtomAddress& o) const {
    return chain_name == o.chain_name && res_id.matches(o.res_id) &&
           atom_name == o.atom_name && altloc == o.altloc;
  }

  std::string str() const {
    return atom_str(chain_name, res_id, atom_name, altloc);
  }
};

} // namespace gemmi

namespace std {
template <> struct hash<gemmi::ResidueId> {
  size_t operator()(const gemmi::ResidueId& r) const {
    size_t seqid_hash = (*r.seqid.num << 7) + (r.seqid.icode | 0x20);
    return seqid_hash ^ hash<string>()(r.segment) ^ hash<string>()(r.name);
  }
};
} // namespace std

#endif
