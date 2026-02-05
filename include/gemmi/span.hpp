//! @file
//! @brief Span and MutableVectorSpan template classes.
//!
//! Span - span of array or std::vector.
//! MutableVectorSpan - span of std::vector with insert() and erase().

// Copyright 2019 Global Phasing Ltd.
//
// Span - span of array or std::vector.
// MutableVectorSpan - span of std::vector with insert() and erase()

#ifndef GEMMI_SPAN_HPP_
#define GEMMI_SPAN_HPP_

#include <algorithm>    // for find_if, find_if_not
#include <vector>
#include <stdexcept>    // for out_of_range
#include <type_traits>  // for remove_cv, conditional, is_const

namespace gemmi {

template<typename Item> struct MutableVectorSpan;

//! @brief Minimalistic span, similar to C++20 std::span.
//! @tparam Item Element type
template<typename Item> struct Span {
  using iterator = Item*;
  using const_iterator = Item const*;
  using element_type = Item;
  using value_type = typename std::remove_cv<Item>::type;

  friend Span<const value_type>;
  friend MutableVectorSpan<value_type>;

  Span() = default;

  //! @brief Constructor from pointer and size.
  //! @param begin Pointer to first element
  //! @param n Number of elements
  Span(iterator begin, std::size_t n) : begin_(begin), size_(n) {}

#if !defined(_MSC_VER) || _MSC_VER-0 >= 1926
  //! @brief Converting constructor (non-const to const).
  //! @tparam T Item type (must be const)
  //! @param o Source span
  template<typename T=Item>
  Span(const Span<value_type>& o,
       typename std::enable_if<std::is_const<T>::value>::type* = 0)
#else
  //! @brief Converting constructor (older MSVC version).
  //! @param o Source span
  Span(const Span<value_type>& o)
#endif
    : begin_(o.begin_), size_(o.size_) {}

  //! @brief Set beginning pointer.
  //! @param begin New beginning pointer
  void set_begin(iterator begin) { begin_ = begin; }

  //! @brief Set size.
  //! @param n New size
  void set_size(std::size_t n) { size_ = n; }

  const_iterator begin() const { return begin_; }  //!< Get const begin iterator
  const_iterator end() const { return begin_ + size_; }  //!< Get const end iterator
  iterator begin() { return begin_; }  //!< Get begin iterator
  iterator end() { return begin_ + size_; }  //!< Get end iterator

  Item& front() { return *begin_; }  //!< Get first element
  const Item& front() const { return *begin_; }  //!< Get first element (const)
  Item& back() { return *(begin_ + size_ - 1); }  //!< Get last element
  const Item& back() const { return *(begin_ + size_ - 1); }  //!< Get last element (const)

  //! @brief Element access operator (const).
  //! @param i Index
  //! @return Reference to element
  const Item& operator[](std::size_t i) const { return *(begin_ + i); }

  //! @brief Element access operator.
  //! @param i Index
  //! @return Reference to element
  Item& operator[](std::size_t i) { return *(begin_ + i); }

  //! @brief Checked element access.
  //! @param i Index
  //! @return Reference to element
  //! @throws std::out_of_range if i >= size()
  Item& at(std::size_t i) {
    if (i >= size())
      throw std::out_of_range("item index ouf of range: #" + std::to_string(i));
    return *(begin_ + i);
  }
  //! @brief Checked element access (const).
  //! @param i Index
  //! @return Reference to element
  //! @throws std::out_of_range if i >= size()
  const Item& at(std::size_t i) const {
    return const_cast<Span*>(this)->at(i);
  }

  std::size_t size() const { return size_; }  //!< Get size
  bool empty() const { return size_ == 0; }  //!< Check if empty
  explicit operator bool() const { return size_ != 0; }  //!< Boolean conversion

  //! @brief Create subspan from iterator range.
  //! @tparam Iter Iterator type
  //! @param first First iterator
  //! @param last Last iterator
  //! @return Subspan
  template<typename Iter> Span<Item> sub(Iter first, Iter last) {
    return Span<Item>(&*first, last - first);
  }

  //! @brief Create subspan matching predicate.
  //! @tparam F Predicate function type
  //! @tparam V Value type (default Item)
  //! @param func Predicate function
  //! @return Subspan of consecutive elements matching func
  template<typename F, typename V=Item> Span<V> subspan(F&& func) {
    iterator group_begin = std::find_if(this->begin(), this->end(), func);
    iterator group_end = std::find_if_not(group_begin, this->end(), func);
    return Span<V>(&*group_begin, group_end - group_begin);
  }
  //! @brief Create const subspan matching predicate.
  //! @tparam F Predicate function type
  //! @param func Predicate function
  //! @return Const subspan
  template<typename F> Span<const value_type> subspan(F&& func) const {
    using V = const value_type;
    return const_cast<Span*>(this)->subspan<F, V>(std::forward<F>(func));
  }

  //! @brief Get children (returns self for iteration).
  //! @return Self reference
  //!
  //! We use children() to iterate over Model, Chain, etc.
  Span& children() { return *this; }

  //! @brief Get children (returns self for iteration, const).
  //! @return Const self reference
  const Span& children() const { return *this; }

private:
  iterator begin_ = nullptr;
  std::size_t size_ = 0;
};

//! @brief Span of std::vector with insert() and erase().
//! @tparam Item Element type
template<typename Item> struct MutableVectorSpan : Span<Item> {
  using vector_type = std::vector<typename Span<Item>::value_type>;  //!< Vector type
  using iterator = typename Span<Item>::iterator;  //!< Iterator type

  MutableVectorSpan() = default;

  //! @brief Constructor from Span and vector pointer.
  //! @param p Span
  //! @param v Vector pointer
  MutableVectorSpan(Span<Item>&& p, vector_type* v)
    : Span<Item>(p), vector_(v) {}

  //! @brief Constructor from vector, begin, and size.
  //! @param v Vector reference
  //! @param begin Begin iterator
  //! @param n Size
  MutableVectorSpan(vector_type& v, iterator begin, std::size_t n)
    : Span<Item>(begin, n), vector_(&v) {}

  //! @brief Create mutable subspan from iterator range.
  //! @tparam Iter Iterator type
  //! @param first First iterator
  //! @param last Last iterator
  //! @return Mutable subspan
  template<typename Iter> MutableVectorSpan<Item> sub(Iter first, Iter last) {
    return {Span<Item>::sub(first, last), vector_};
  }

  //! @brief Create mutable subspan matching predicate.
  //! @tparam F Predicate function type
  //! @param func Predicate function
  //! @return Mutable subspan
  template<typename F> MutableVectorSpan<Item> subspan(F&& func) {
    return {Span<Item>::subspan(std::forward<F>(func)), vector_};
  }

  //! @brief Create const mutable subspan matching predicate.
  //! @tparam F Predicate function type
  //! @param func Predicate function
  //! @return Const mutable subspan
  template<typename F> MutableVectorSpan<const Item> subspan(F&& func) const {
    return {Span<const Item>::subspan(std::forward<F>(func)), vector_};
  }

  //! @brief Insert item at position.
  //! @param pos Iterator position
  //! @param item Item to insert (moved)
  //! @return Iterator to inserted item
  iterator insert(iterator pos, Item&& item) {
    auto offset = this->begin_ - this->vector_->data();
    auto iter = vector_->begin() + (pos - this->vector_->data());
    auto ret = vector_->insert(iter, std::move(item));
    this->begin_ = vector_->data() + offset;
    ++this->size_;
    return &*ret;
  }

  //! @brief Erase item at position.
  //! @param pos Iterator position
  void erase(iterator pos) {
    vector_->erase(vector_->begin() + (pos - vector_->data()));
    --this->size_;
  }

  //! @brief Check if span starts at vector beginning.
  //! @return True if at beginning
  bool is_beginning() const { return this->begin() == vector_->data(); }

  //! @brief Check if span ends at vector end.
  //! @return True if at end
  bool is_ending() const { return this->end() == vector_->data() + vector_->size(); }

private:
  vector_type* vector_ = nullptr;  //!< Vector pointer for insert() and erase()
};

} // namespace gemmi
#endif
