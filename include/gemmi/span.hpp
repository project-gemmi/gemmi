// Copyright 2019 Global Phasing Ltd.
//
// VectorSpan - span of std::vector elements.

#ifndef GEMMI_SPAN_HPP_
#define GEMMI_SPAN_HPP_

#include <vector>
#include <stdexcept>    // for out_of_range
#include <type_traits>  // for remove_cv

namespace gemmi {

// It is somewhat similar to C++20 std::span, but works only with std::vector,
// and implements insert() and erase().
template<typename Item> struct VectorSpan {
  using iterator = Item*;
  using const_iterator = Item const*;
  using value_type = typename std::remove_cv<Item>::type;

  VectorSpan() = default;
  VectorSpan(std::vector<Item>& v, iterator begin, std::size_t n)
    : begin_(begin), size_(n), vector_(&v) {}
  VectorSpan(VectorSpan& base, iterator begin, std::size_t n)
    : begin_(begin), size_(n), vector_(base.vector_) {}
  VectorSpan(std::vector<Item>& v)
    : begin_(v.data()), size_(v.size()), vector_(&v) {}

  const_iterator begin() const { return begin_; }
  const_iterator end() const { return begin_ + size_; }
  iterator begin() { return begin_; }
  iterator end() { return begin_ + size_; }

  Item& front() { return *begin_; }
  const Item& front() const { return *begin_; }
  Item& back() { return *(begin_ + size_ - 1); }
  const Item& back() const { return *(begin_ + size_ - 1); }

  const Item& operator[](std::size_t i) const { return *(begin_ + i); }
  Item& operator[](std::size_t i) { return *(begin_ + i); }

  Item& at(std::size_t i) {
    if (i >= size())
      throw std::out_of_range("item index ouf of range: #" + std::to_string(i));
    return *(begin_ + i);
  }
  const Item& at(std::size_t i) const {
    return const_cast<VectorSpan*>(this)->at(i);
  }

  std::size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  explicit operator bool() const { return size_ != 0; }

  iterator insert(iterator pos, Item&& item) {
    auto offset = begin_ - vector_->data();
    auto iter = vector_->begin() + (pos - vector_->data());
    auto ret = vector_->insert(iter, std::move(item));
    begin_ = vector_->data() + offset;
    ++size_;
    return &*ret;
  }

  void erase(iterator pos) {
    vector_->erase(vector_->begin() + (pos - vector_->data()));
    --size_;
  }

private:
  iterator begin_;
  std::size_t size_ = 0;
  std::vector<Item>* vector_ = nullptr;  // for insert() and erase()
};

} // namespace gemmi
#endif
