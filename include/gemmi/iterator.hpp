// Copyright 2018 Global Phasing Ltd.
//
// Iterators. Currently each of them is a BidirectionalIterator.

#ifndef GEMMI_ITERBASE_HPP_
#define GEMMI_ITERBASE_HPP_
#include <iterator>  // for bidirectional_iterator_tag
#include <utility>  // for forward
#include <vector>

namespace gemmi {

// implements concept BidirectionalIterator
template <typename Policy>
struct BidirIterator : Policy {
  typedef typename Policy::value_type value_type;
  typedef std::ptrdiff_t difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef std::bidirectional_iterator_tag iterator_category;

  BidirIterator() = default;
  BidirIterator(Policy&& p) : Policy(p) {}

  BidirIterator &operator++() { Policy::increment(); return *this; }
  BidirIterator operator++(int) { BidirIterator x = *this; ++*this; return x; }
  BidirIterator &operator--() { Policy::decrement(); return *this; }
  BidirIterator operator--(int) { BidirIterator x = *this; --*this; return x; }
  bool operator==(const BidirIterator &o) const { return Policy::equal(o); }
  bool operator!=(const BidirIterator &o) const { return !Policy::equal(o); }
  reference operator*() { return Policy::dereference(); }
  pointer operator->() { return &Policy::dereference(); }
  using const_variant = BidirIterator<typename Policy::const_policy>;
  operator const_variant() const {
    return const_variant(static_cast<const Policy&>(*this));
  }
};

template<typename Value>
class StrideIterPolicy {
public:
  typedef Value value_type;
  StrideIterPolicy() : cur_(nullptr), offset_(0), stride_(0) {}
  StrideIterPolicy(Value* ptr, std::size_t offset, unsigned stride)
    : cur_(ptr), offset_(offset), stride_(stride) {}
  void increment() { cur_ += stride_; }
  void decrement() { cur_ -= stride_; }
  bool equal(const StrideIterPolicy& o) const { return cur_ == o.cur_; }
  Value& dereference() { return cur_[offset_]; }
  using const_policy = StrideIterPolicy<Value const>;
  operator const_policy() const { return const_policy(cur_, offset_, stride_); }
private:
  Value* cur_;
  unsigned offset_;
  unsigned stride_;
};
template<typename Value>
using StrideIter = BidirIterator<StrideIterPolicy<Value>>;


template<typename Redirect, typename Value>
class IndirectIterPolicy {
public:
  typedef Value value_type;
  IndirectIterPolicy() : redir_(nullptr) {}
  IndirectIterPolicy(Redirect* redir, std::vector<int>::const_iterator cur)
    : redir_(redir), cur_(cur) {}
  void increment() { ++cur_; }
  void decrement() { --cur_; }
  bool equal(const IndirectIterPolicy& o) const { return cur_ == o.cur_; }
  Value& dereference() { return redir_->value_at(*cur_); }
  using const_policy = IndirectIterPolicy<Redirect const, Value const>;
  operator const_policy() const { return const_policy(redir_, cur_); }
  // TODO: what should be done with absent optional tags (*cur_ < 0)?
private:
  Redirect* redir_;
  std::vector<int>::const_iterator cur_; // points into positions
};
template<typename Redirect, typename Value>
using IndirectIter = BidirIterator<IndirectIterPolicy<Redirect, Value>>;


template<typename Value>
class UniqIterPolicy {
public:
  typedef Value value_type;
  UniqIterPolicy() : vec_(nullptr), pos_(0) {}
  UniqIterPolicy(std::vector<Value>* vec, std::size_t pos)
    : vec_(vec), pos_(pos) {}
  void increment() { ++pos_; } // TODO
  void decrement() { --pos_; }
  bool equal(const UniqIterPolicy& o) const { return pos_ == o.pos_; }
  Value& dereference() { return vec_[pos_]; }
  using const_policy = UniqIterPolicy<Value const>;
  operator const_policy() const { return const_policy(vec_, pos_); }
private:
  std::vector<Value>* vec_;
  std::size_t pos_;
};
template<typename Value>
using UniqIter = BidirIterator<UniqIterPolicy<Value>>;

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
