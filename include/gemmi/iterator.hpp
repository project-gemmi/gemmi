// Copyright 2018 Global Phasing Ltd.
//
// Iterators. Currently each of them is a BidirectionalIterator.

#ifndef GEMMI_ITERATOR_HPP_
#define GEMMI_ITERATOR_HPP_
#include <iterator>     // for bidirectional_iterator_tag
#include <type_traits>  // for remove_cv
#include <vector>

#ifdef  __INTEL_COMPILER
// warning #597: "X<T>::operator X<T>() const" will not be called for implicit
// or explicit conversions. That warning is triggered when templates
// StrideIter, IndirectIter and others are expanded with const Value.
# pragma warning disable 597
#endif

namespace gemmi {

// implements concept BidirectionalIterator
template <typename Policy>
struct BidirIterator : Policy {
  typedef typename std::remove_cv<typename Policy::value_type>::type value_type;
  typedef std::ptrdiff_t difference_type;
  typedef typename Policy::value_type* pointer;
  typedef typename Policy::value_type& reference;
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
  std::size_t offset_;
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


template<typename Vector, typename Value>
class UniqIterPolicy {
public:
  typedef Value value_type;
  UniqIterPolicy() : vec_(nullptr), pos_(0) {}
  UniqIterPolicy(Vector* vec, std::size_t pos) : vec_(vec), pos_(pos) {}
  void increment() {
    std::size_t old = pos_++;
    while (pos_ != vec_->size() && (*vec_)[old].same_group((*vec_)[pos_]))
      ++pos_;
  }
  void decrement() {
    --pos_;
    while (pos_ != 0 && (*vec_)[pos_ - 1].same_group((*vec_)[pos_]))
      --pos_;
  }
  bool equal(const UniqIterPolicy& o) const { return pos_ == o.pos_; }
  Value& dereference() { return (*vec_)[pos_]; }
  using const_policy = UniqIterPolicy<Vector const, Value const>;
  operator const_policy() const { return const_policy(vec_, pos_); }
private:
  Vector* vec_;
  std::size_t pos_;
};
template<typename Vector, typename Value>
using UniqIter = BidirIterator<UniqIterPolicy<Vector, Value>>;

template<typename Value, typename Vector=std::vector<Value>>
struct UniqProxy {
  Vector& vec;
  using iterator = UniqIter<Vector, Value>;
  iterator begin() { return {{&vec, 0}}; }
  iterator end() { return {{&vec, vec.size()}}; }
};
template<typename Value, typename Vector=std::vector<Value>>
struct ConstUniqProxy {
  const Vector& vec;
  using iterator = UniqIter<const Vector, const Value>;
  iterator begin() const { return {{&vec, 0}}; }
  iterator end() const { return {{&vec, vec.size()}}; }
};


template<typename Filter, typename Vector, typename Value>
class FilterIterPolicy {
public:
  typedef Value value_type;
  FilterIterPolicy() : vec_(nullptr), pos_(0) {}
  FilterIterPolicy(const Filter* filter, Vector* vec, std::size_t pos)
      : filter_(filter), vec_(vec), pos_(pos) {
    while (pos_ != vec_->size() && !matches(pos_)) ++pos_;
  }
  bool matches(std::size_t p) const { return filter_->matches((*vec_)[p]); }
  void increment() { while (pos_ != vec_->size() && !matches(++pos_)) {} }
  void decrement() { while (pos_ != 0 && !matches(--pos_)) {} }
  bool equal(const FilterIterPolicy& o) const { return pos_ == o.pos_; }
  Value& dereference() { return (*vec_)[pos_]; }
  using const_policy = FilterIterPolicy<Filter, Vector const, Value const>;
  operator const_policy() const { return const_policy(vec_, pos_); }
private:
  const Filter* filter_;
  Vector* vec_;
  std::size_t pos_;
};
template<typename Filter, typename Vector, typename Value>
using FilterIter = BidirIterator<FilterIterPolicy<Filter, Vector, Value>>;

template<typename Filter, typename Value>
struct FilterProxy {
  const Filter& filter;
  std::vector<Value>& vec;
  using iterator = FilterIter<Filter, std::vector<Value>, Value>;
  iterator begin() { return {{&filter, &vec, 0}}; }
  iterator end() { return {{&filter, &vec, vec.size()}}; }
};

template<typename Filter, typename Value>
struct ConstFilterProxy {
  const Filter& filter;
  const std::vector<Value>& vec;
  using iterator = FilterIter<Filter, const std::vector<Value>, const Value>;
  iterator begin() const { return {{&filter, &vec, 0}}; }
  iterator end() const { return {{&filter, &vec, vec.size()}}; }
};

} // namespace gemmi
#endif
