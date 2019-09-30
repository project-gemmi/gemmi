// Copyright 2017 Global Phasing Ltd.
//
// SeqId -- residue number and insertion code together

#ifndef GEMMI_SEQID_HPP_
#define GEMMI_SEQID_HPP_

#include <cstdlib>    // for strtol
#include <stdexcept>  // for invalid_argument
#include <string>

namespace gemmi {

// Optional int value. N is a special value that means not-set.
template<int N> struct OptionalInt {
  enum { None=N };
  int value = None;

  OptionalInt() = default;
  OptionalInt(int n) : value(n) {}
  bool has_value() const { return value != None; }
  std::string str(char null='?') const {
    return has_value() ? std::to_string(value) : std::string(1, null);
  }
  OptionalInt& operator=(int n) { value = n; return *this; }
  bool operator==(const OptionalInt& o) const { return value == o.value; }
  bool operator!=(const OptionalInt& o) const { return value != o.value; }
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
};

struct SeqId {
  using OptionalNum = OptionalInt<-999>;

  OptionalNum num; // sequence number
  char icode = ' ';  // insertion code

  SeqId() = default;
  SeqId(int num_, char icode_) { num = num_; icode = icode_; }
  explicit SeqId(const std::string& str) {
    char* endptr;
    num = std::strtol(str.c_str(), &endptr, 10);
    if (endptr == str.c_str() || (*endptr != '\0' && endptr[1] != '\0'))
      throw std::invalid_argument("Not a seqid: " + str);
    icode = (*endptr | 0x20);
  }

  bool operator==(const SeqId& o) const {
    return num == o.num && (icode | 0x20) == (o.icode | 0x20);
  }
  bool operator!=(const SeqId& o) const { return !operator==(o); }

  char has_icode() const { return icode != ' '; }

  std::string str() const {
    std::string r = num.str();
    if (icode != ' ')
      r += icode;
    return r;
  }
};

} // namespace gemmi
#endif
