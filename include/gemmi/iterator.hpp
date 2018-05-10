// Copyright 2018 Global Phasing Ltd.
//
// Iterators. Currently each of them is a BidirectionalIterator.

#ifndef GEMMI_ITERBASE_HPP_
#define GEMMI_ITERBASE_HPP_
#include <iterator>  // for bidirectional_iterator_tag
#include <utility>  // for forward
#include <vector>

namespace gemmi {
namespace cif {

// implements concept BidirectionalIterator
template <typename Value, typename Policy>
struct BidirIterator : Policy {
  typedef Value value_type;
  typedef std::ptrdiff_t difference_type;
  typedef Value* pointer;
  typedef Value& reference;
  typedef std::bidirectional_iterator_tag iterator_category;

  BidirIterator() = default;
  template<typename... Args>
  BidirIterator(Args&&... args) : Policy(std::forward<Args>(args)...) {}

  BidirIterator &operator++() { Policy::increment(); return *this; }
  BidirIterator operator++(int) { auto x = *this; ++*this; return x; }
  BidirIterator &operator--() { Policy::decrement(); return *this; }
  BidirIterator operator--(int) { auto x = *this; --*this; return x; }
  bool operator==(const BidirIterator &o) const { return Policy::equal(o); }
  bool operator!=(const BidirIterator &o) const { return !Policy::equal(o); }
  Value& operator*() { return Policy::dereference(); }
  Value* operator->() { return &Policy::dereference(); }
};

template<typename Value>
class StrideIterPolicy {
public:
  StrideIterPolicy() : cur_(nullptr), offset_(0), stride_(0) {}
  StrideIterPolicy(Value* ptr, std::size_t offset, unsigned stride)
    : cur_(ptr), offset_(offset), stride_(stride) {}
  void increment() { cur_ += stride_; }
  void decrement() { cur_ -= stride_; }
  bool equal(const StrideIterPolicy& o) const { return cur_ == o.cur_; }
  Value& dereference() { return cur_[offset_]; }
  operator StrideIterPolicy<Value const>() const {
    return StrideIterPolicy<Value const>(cur_, offset_, stride_);
  }
private:
  Value* cur_;
  unsigned offset_;
  unsigned stride_;
};
template<typename Value>
using StrideIter = BidirIterator<Value, StrideIterPolicy<Value>>;


template<typename Redirect, typename Value>
class IndirectIterPolicy {
public:
  IndirectIterPolicy() : redir_(nullptr) {}
  IndirectIterPolicy(Redirect* redir, std::vector<int>::const_iterator cur)
    : redir_(redir), cur_(cur) {}
  void increment() { ++cur_; }
  void decrement() { --cur_; }
  bool equal(const IndirectIterPolicy& o) const { return cur_ == o.cur_; }
  Value& dereference() { return redir_->value_at(*cur_); }
  operator IndirectIterPolicy<Redirect const, Value const>() const {
    return IndirectIterPolicy<Redirect const, Value const>(redir_, cur_);
  }
  // TODO: what should be done with absent optional tags (*cur_ < 0)?
private:
  Redirect* redir_;
  std::vector<int>::const_iterator cur_; // points into positions
};
template<typename Redirect, typename Val>
using IndirectIter = BidirIterator<Val, IndirectIterPolicy<Redirect, Val>>;

} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
