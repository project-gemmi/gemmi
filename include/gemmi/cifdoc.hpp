// Copyright 2017 Global Phasing Ltd.
//
// struct Document that represents the CIF file (but could be also
// read from, for example, JSON file (CIF-JSON or mmJSON).

#ifndef GEMMI_CIFDOC_HPP_
#define GEMMI_CIFDOC_HPP_
#include "util.hpp"  // for starts_with, to_lower
#include <algorithm> // for move, find_if, all_of, min
#include <array>
#include <cassert>
#include <cstring>   // for memchr
#include <initializer_list>
#include <iosfwd>    // for size_t, ptrdiff_t
#include <iterator>  // for bidirectional_iterator_tag
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace gemmi {
namespace cif {
using std::size_t;

// base for a BidirectionalIterator (std::iterator is deprecated in C++17)
template <typename Value> struct IterBase {
  typedef Value value_type;
  typedef std::ptrdiff_t  difference_type;
  typedef Value* pointer;
  typedef Value& reference;
  typedef std::bidirectional_iterator_tag iterator_category;
};

enum class ItemType : unsigned char {
  Value,
  Loop,
  Frame,
  Comment,
  Erased,
};

enum class ValueType : unsigned char {
  NotSet,
  Char, // Line/Text?
  Numb, // Int/Float?
  Dot,
  QuestionMark,
};

inline bool is_null(const std::string& value) {
  return value == "?" || value == ".";
}

inline std::string value_type_to_str(ValueType v) {
  switch (v) {
    case ValueType::NotSet: return "n/a";
    case ValueType::Char: return "char";
    case ValueType::Numb: return "numb";
    case ValueType::Dot: return "'.'";
    case ValueType::QuestionMark: return "'?'";
  }
  return "";
}

inline std::string as_string(const std::string& value) {
  if (value.empty() || value == "?" || value == ".")
    return "";
  if (value[0] == '"' || value[0] == '\'')
    return std::string(value.begin() + 1, value.end() - 1);
  if (value[0] == ';' && value.size() > 2 && *(value.end() - 2) == '\n')
    return std::string(value.begin() + 1, value.end() - 2);
  return value;
}

inline std::string as_string(const std::string* value) {
  return value ? as_string(*value) : std::string();
}

using Pair = std::array<std::string, 2>;

// used only as arguments when creating Item
struct LoopArg {};
struct FrameArg { std::string str; };
struct CommentArg { std::string str; };

struct LoopTag {
  ValueType valtype = ValueType::NotSet;
  int line_number = -1;
  std::string tag;
  explicit LoopTag(std::string&& t) : tag{t} {}
};

struct Loop {
  std::vector<LoopTag> tags;
  std::vector<std::string> values;

  // search and access
  int find_tag(const std::string& tag) const {
    auto f = std::find_if(tags.begin(), tags.end(),
                          [&tag](const LoopTag& lt) { return lt.tag == tag; });
    return f == tags.end() ? -1 : f - tags.begin();
  }
  int must_find_tag(const std::string& tag) const {
    auto f = std::find_if(tags.begin(), tags.end(),
                          [&tag](const LoopTag& lt) { return lt.tag == tag; });
    if (f == tags.end())
      throw std::runtime_error("required tag not found: " + tag);
    return f - tags.begin();
  }
  size_t width() const { return tags.size(); }
  size_t length() const { return values.size() / tags.size(); }
  const std::string& val(size_t row, size_t col) const {
    return values[row * tags.size() + col];
  }
  void clear() { tags.clear(); values.clear(); }
  void append_row(std::initializer_list<std::string>&& vv) {
    assert(vv.size() == tags.size());
    for (const std::string& v : vv)
      values.emplace_back(v);
  }
};


template<typename T>
class StrideIter : public IterBase<T> {
public:
  StrideIter() : cur_(nullptr), offset_(0), stride_(0) {}
  StrideIter(T* ptr, size_t offset, unsigned stride)
    : cur_(ptr), offset_(offset), stride_(stride) {}
  StrideIter& operator++() { cur_ += stride_; return *this; }
  StrideIter& operator--() { cur_ -= stride_; return *this; }
  StrideIter operator++(int) { auto t = *this; cur_ += stride_; return t; }
  StrideIter operator--(int) { auto t = *this; cur_ -= stride_; return t; }
  const T& operator*() const { return cur_[offset_]; }
  const T* operator->() const { return cur_ + offset_; }
  T& operator*() { return cur_[offset_]; }
  T* operator->() { return cur_ + offset_; }
  bool operator==(const StrideIter& other) const { return cur_ == other.cur_; }
  bool operator!=(const StrideIter& other) const { return cur_ != other.cur_; }
private:
  T* cur_;
  unsigned offset_;
  unsigned stride_;
};

struct Item;
struct Block;

// Accessor to a specific loop column, or to a single value from a Pair.
struct Column {
  const Item* it;
  size_t col;  // for loop it is a column index in it->loop

  using iterator = StrideIter<std::string>;
  using const_iterator = StrideIter<const std::string>;
  const_iterator begin() const;
  const_iterator end() const;
  //iterator begin();
  //iterator end();

  const Loop* get_loop() const;
  const std::string* get_tag() const;
  int length() const {
    if (const Loop* loop = get_loop())
      return loop->length();
    return it ? 1 : 0;
  }
  const std::string& operator[](int n) const;
  const std::string& at(int n) const {
    if (n < 0)
      n += length();
    if (n < 0 || n >= length())
      throw std::out_of_range("Cannot access element " + std::to_string(n) +
          " in Column with length " + std::to_string(length()));
    return operator[](n);
  }
  std::string str(int n) const { return as_string(at(n)); }
};

// Some values can be given either in loop or as tag-value pairs.
// The latter case is equivalent to a loop with a single row.
// We optimized for loops, and in case of tag-values we copy the values
// into the `values` vector.
struct Table {
  const Loop* loop;
  const Block& blo;
  std::vector<int> positions;

  struct Row {
    const Table& tab;
    int row_index;

    const std::string& direct_at(int pos) const;
    const std::string& at(int n) const {
      if (n < 0)
        n += size();
      return direct_at(tab.positions.at(n));
    }
    const std::string& operator[](int n) const {
      return direct_at(tab.positions[n]);
    }
    const std::string* ptr_at(int n) const {
      if (n < 0)
        n += size();
      int pos = tab.positions.at(n);
      return pos >= 0 ? &direct_at(pos) : nullptr;
    }
    bool has(int n) const { return tab.has_column(n); }
    size_t size() const { return tab.width(); }
    std::string str(int n) const { return as_string(at(n)); }
    struct iterator : IterBase<std::string> {
      iterator() : parent(nullptr) {}
      iterator(const Row* p, std::vector<int>::const_iterator it)
        : parent(p), cur(it) {}
      const Row* parent;
      std::vector<int>::const_iterator cur; // points into positions
      iterator& operator++() { ++cur; return *this; }
      iterator& operator--() { --cur; return *this; }
      iterator operator++(int) { auto t = *this; ++cur; return t; }
      iterator operator--(int) { auto t = *this; --cur; return t; }
      const std::string& operator*() const { return parent->direct_at(*cur); }
      const std::string* operator->() const { return &(operator*()); }
      bool operator==(const iterator& other) const { return cur == other.cur; }
      bool operator!=(const iterator& other) const { return cur != other.cur; }
      // TODO: what should be done with absent optional tags?
    };
    using const_iterator = iterator;
    const_iterator begin() const {
      return const_iterator(this, tab.positions.begin());
    }
    const_iterator end() const {
      return const_iterator(this, tab.positions.end());
    }
    iterator begin() { return iterator(this, tab.positions.begin()); }
    iterator end() { return iterator(this, tab.positions.end()); }
  };

  bool ok() const { return !positions.empty(); }
  size_t width() const { return positions.size(); }
  size_t length() const {
    return loop ? loop->length() : (positions.empty() ? 0 : 1);
  }
  Row tags() const { return Row{*this, -1}; }
  bool has_column(int n) const { return positions.at(n) >= 0; }

  Row operator[](int n) const { return Row{*this, n}; }

  Row at(int n) const {
    if (n < 0)
      n += length();
    if (n < 0 || static_cast<size_t>(n) >= length())
      throw std::out_of_range("No row with index " + std::to_string(n));
    return (*this)[n];
  }

  Row one() const {
    if (length() != 1)
      throw std::runtime_error("Expected one value, found " +
                                std::to_string(length()));
    return (*this)[0];
  }

  Row find_row(const std::string& s) const;
  // TODO: get_column(n), find_column(name)

  // It is not a proper input iterator, but just enough for using range-for.
  struct iterator {
    const Table& parent;
    size_t index;
    void operator++() { index++; }
    bool operator==(const iterator& o) const { return index == o.index; }
    bool operator!=(const iterator& o) const { return index != o.index; }
    Row operator*() const { return parent[index]; }
    const std::string& get(int n) const { return parent[index].at(n); }
  };
  iterator begin() const { return iterator{*this, 0}; }
  iterator end() const { return iterator{*this, length()}; }
};

struct Block {
  std::string name;
  std::vector<Item> items;

  explicit Block(const std::string& name_) : name(name_) {}
  Block() {}

  // access functions
  const Item* find_pair_item(const std::string& tag) const;
  const Pair* find_pair(const std::string& tag) const;
  const std::string* find_value(const std::string& tag) const {
    const Pair* pair = find_pair(tag);
    return pair ? &(*pair)[1] : nullptr;
  }
  Column find_loop(const std::string& tag) const;
  Column find_values(const std::string& tag) const;
  Table find(const std::string& prefix,
             const std::vector<std::string>& tags) const;
  Table find(const std::vector<std::string>& tags) const {
    return find({}, tags);
  }

  // modifying functions
  void set_pair(const std::string& tag, std::string v);
  bool delete_loop(const std::string& tag);
  // These functions delete all keys/loops that start with the prefix.
  // For mmCIF the prefix should normally end with dot.
  Loop& clear_or_add_loop(const std::string& prefix,
                          const std::initializer_list<const char*>& tags);
  void delete_category(const std::string& prefix);

  // mmCIF specific functions
  std::vector<std::string> get_mmcif_category_names() const;
  Table find_mmcif_category(std::string cat);
};

struct Item {
  ItemType type;
  ValueType valtype = ValueType::NotSet; // for Pair only
  int line_number = -1;
  union {
    Pair pair;
    Loop loop;
    Block frame;
  };

  explicit Item(LoopArg)
    : type{ItemType::Loop}, loop{} {}
  explicit Item(std::string&& t)
    : type{ItemType::Value}, pair{std::move(t), std::string()} {}
  Item(const std::string& t, const std::string& v)
    : type{ItemType::Value}, pair{t, v} {}
  explicit Item(FrameArg&& frame_arg)
    : type{ItemType::Frame}, frame(frame_arg.str) {}
  explicit Item(CommentArg&& comment)
    : type{ItemType::Comment}, pair{std::string(), std::move(comment.str)} {}

  Item(Item&& o) noexcept
      : type(o.type), valtype(o.valtype), line_number(o.line_number) {
    move_value(std::move(o));
  }
  Item(const Item&& o)
      : type(o.type), valtype(o.valtype), line_number(o.line_number) {
    copy_value(o);
  }
  Item(const Item& o)
      : type(o.type), valtype(o.valtype), line_number(o.line_number) {
    copy_value(o);
  }
  Item(Item& o) : Item(static_cast<const Item&>(o)) {}

  ~Item() {
    switch (type) {
      case ItemType::Value: pair.~Pair(); break;
      case ItemType::Loop: loop.~Loop(); break;
      case ItemType::Frame: frame.~Block(); break;
      case ItemType::Comment: pair.~Pair(); break;
      case ItemType::Erased: break;
    }
  }

  void erase() {
    this->~Item();
    type = ItemType::Erased;
  }

  bool has_prefix(const std::string& prefix) const {
    return (type == ItemType::Value && gemmi::starts_with(pair[0], prefix)) ||
           (type == ItemType::Loop && !loop.tags.empty() &&
            gemmi::starts_with(loop.tags[0].tag, prefix));
  }

private:
  void copy_value(const Item& o) {
    if (o.type == ItemType::Value || o.type == ItemType::Comment)
      new (&pair) Pair(o.pair);
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(o.loop);
    else if (o.type == ItemType::Frame)
      new (&frame) Block(o.frame);
  }

  void move_value(Item&& o) {
    if (o.type == ItemType::Value || o.type == ItemType::Comment)
      new (&pair) Pair(std::move(o.pair));
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(std::move(o.loop));
    else if (o.type == ItemType::Frame)
      new (&frame) Block(std::move(o.frame));
  }
};

inline const std::string* Column::get_tag() const {
  if (!it)
    return nullptr;
  if (const Loop* loop = get_loop())
    return &loop->tags.at(col).tag;
  return &it->pair[0];
}

inline const Loop* Column::get_loop() const {
  return it && it->type == ItemType::Loop ? &it->loop : nullptr;
}
inline Column::const_iterator Column::begin() const {
  if (const Loop* loop = get_loop())
    return const_iterator(loop->values.data(), col, loop->width());
  if (it && it->type == ItemType::Value)
    return const_iterator(&it->pair[1], 0, 1);
  return const_iterator();
}

inline Column::const_iterator Column::end() const {
  if (const Loop* loop = get_loop())
    return const_iterator(loop->values.data() + loop->values.size(),
                          col, loop->width());
  if (it && it->type == ItemType::Value)
    return const_iterator(&it->pair[1] + 1, 0, 1);
  return const_iterator();
}

inline const std::string& Column::operator[](int n) const {
  if (const Loop* loop = get_loop())
    return loop->values[n * loop->width() + col];
  return it->pair[1];
}

inline const std::string& Table::Row::direct_at(int pos) const {
  if (pos == -1)
    throw std::out_of_range("Cannot access missing optional tag.");
  if (row_index == -1) { // tags
    if (tab.loop)
      return tab.loop->tags.at(pos).tag;
    return tab.blo.items[pos].pair[0];
  }
  if (tab.loop)
    return tab.loop->values.at(tab.loop->width() * row_index + pos);
  return tab.blo.items[pos].pair[1];
}

inline Table::Row Table::find_row(const std::string& s) const {
  int pos = positions.at(0);
  if (loop) {
    for (size_t i = 0; i < loop->values.size(); i += loop->width())
      if (as_string(loop->values[i + pos]) == s)
        return Row{*this, static_cast<int>(i / loop->width())};
  } else if (as_string(blo.items[pos].pair[1]) == s) {
    return Row{*this, 0};
  }
  throw std::runtime_error("Not found in the first column: " + s);
}

inline const Item* Block::find_pair_item(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.pair[0] == tag)
      return &i;
  return nullptr;
}

inline const Pair* Block::find_pair(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.pair[0] == tag)
      return &i.pair;
  return nullptr;
}

inline void Block::set_pair(const std::string& tag, std::string v) {
  if (tag[0] != '_')
    throw std::runtime_error("Tag should start with '_', got: " + tag);
  for (Item& i : items) {
    if (i.type == ItemType::Value && i.pair[0] == tag) {
      i.pair[1] = v;
      return;
    }
    if (i.type == ItemType::Loop && i.loop.find_tag(tag) != -1) {
      i.loop.~Loop();
      i.type = ItemType::Value;
      i.pair = { tag, v };
      return;
    }
  }
  items.emplace_back(tag, v);
}

inline Column Block::find_loop(const std::string& tag) const {
  Column c = find_values(tag);
  return c.it && c.it->type == ItemType::Loop ? c : Column{nullptr, 0};
}

inline Column Block::find_values(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Loop) {
      int pos = i.loop.find_tag(tag);
      if (pos != -1)
        return Column{&i, static_cast<size_t>(pos)};
    } else if (i.type == ItemType::Value) {
      if (i.pair[0] == tag)
        return Column{&i, 0};
    }
  return Column{nullptr, 0};
}

inline bool Block::delete_loop(const std::string& tag) {
  for (Item& i : items)
    if (i.type == ItemType::Loop && i.loop.find_tag(tag) != -1) {
      i.erase();
      return true;
    }
  return false;
}

inline void Block::delete_category(const std::string& prefix) {
  for (Item& i : items)
    if (i.has_prefix(prefix))
      i.erase();
}

inline std::vector<std::string> Block::get_mmcif_category_names() const {
  std::vector<std::string> cats;
  for (const Item& item : items) {
    const std::string* tag = nullptr;
    if (item.type == ItemType::Value)
      tag = &item.pair[0];
    else if (item.type == ItemType::Loop && !item.loop.tags.empty())
      tag = &item.loop.tags[0].tag;
    if (tag)
      for (auto j = cats.rbegin(); j != cats.rend(); ++j)
        if (gemmi::starts_with(*tag, *j)) {
          tag = nullptr;
          break;
        }
    if (tag) {
      size_t dot = tag->find('.');
      if (dot != std::string::npos)
        cats.emplace_back(*tag, 0, dot + 1);
    }
  }
  return cats;
}

inline Loop& Block::clear_or_add_loop(const std::string& prefix,
                              const std::initializer_list<const char*>& tags) {
  for (Item& i : items)
    if (i.type == ItemType::Loop && i.has_prefix(prefix)) {
      i.loop.clear();
      return i.loop;
    }
  delete_category(prefix);
  items.emplace_back(LoopArg{});
  Loop& loop = items.back().loop;
  loop.tags.reserve(tags.size());
  for (const char* tag : tags)
    loop.tags.emplace_back(prefix + tag);
  return loop;
}

inline Table Block::find(const std::string& prefix,
                         const std::vector<std::string>& tags) const {
  const Loop* loop = nullptr;
  if (!tags.empty())
    if (const Item* item = find_loop(prefix + tags[0]).it)
      loop = &item->loop;

  std::vector<int> indices;
  indices.reserve(tags.size());
  if (loop) {
    for (const std::string& tag : tags) {
      int idx = loop->find_tag(prefix + (tag[0] != '?' ? tag : tag.substr(1)));
      if (idx == -1 && tag[0] != '?') {
        loop = nullptr;
        indices.clear();
        break;
      }
      indices.push_back(idx);
    }
  } else {
    for (const std::string& tag : tags) {
      // TODO: ?optional_tag
      if (const Item* p = find_pair_item(prefix + tag)) {
        indices.push_back(p - items.data());
      } else {
        indices.clear();
        break;
      }
    }
  }
  return Table{loop, *this, indices};
}

inline Table Block::find_mmcif_category(std::string cat) {
  if (cat[0] != '_')
    throw std::runtime_error("Category should start with '_', got: " + cat);
  if (*(cat.end() - 1) != '.')
    cat += '.';
  std::vector<int> indices;
  for (const Item& i : items)
    if (i.has_prefix(cat)) {
      if (i.type == ItemType::Loop) {
        std::vector<int> indices(i.loop.tags.size());
        for (size_t j = 0; j != indices.size(); ++j) {
          indices[j] = j;
          const std::string& tag = i.loop.tags[j].tag;
          if (!starts_with(tag, cat))
            throw std::runtime_error("Tag " + tag + " in loop with " + cat);
        }
        return Table{&i.loop, *this, indices};
      } else {
        indices.push_back(&i - items.data());
      }
    }
  return Table{nullptr, *this, indices};
}


struct Document {
  std::string source;
  std::vector<Block> blocks;

  // implementation detail: items of the currently parsed block or frame
  std::vector<Item>* items_ = nullptr;

  void clear() noexcept {
    source.clear();
    blocks.clear();
    items_ = nullptr;
  }

  // returns blocks[0] if the document has exactly one block (like mmCIF)
  const Block& sole_block() const {
    if (blocks.size() > 1)
      throw std::runtime_error(std::to_string(blocks.size()) + " data blocks,"
                               " a single block was expected");
    return blocks.at(0);
  }

  const Block* find_block(const std::string& name) const {
    for (const Block& b : blocks)
      if (b.name == name)
        return &b;
    return nullptr;
  }
};


[[noreturn]]
inline void cif_fail(const Document& d, const Block& b, const Item& item,
                     const std::string& s) {
  throw std::runtime_error(d.source + ":" + std::to_string(item.line_number) +
                           " in data_" + b.name + ": " + s);
}

// Throw an error if any block name, frame name or tag is duplicated.
inline void check_duplicates(const Document& d) {
  // check for duplicate block names (except empty "" which is global_)
  std::unordered_set<std::string> names;
  for (const Block& block : d.blocks) {
    bool ok = names.insert(gemmi::to_lower(block.name)).second;
    if (!ok && !block.name.empty())
      throw std::runtime_error("duplicate block name: " + block.name);
  }
  // check for dups inside each block
  std::unordered_set<std::string> frame_names;
  for (const Block& block : d.blocks) {
    names.clear();
    frame_names.clear();
    for (const Item& item : block.items) {
      if (item.type == ItemType::Value) {
        bool ok = names.insert(gemmi::to_lower(item.pair[0])).second;
        if (!ok)
          cif_fail(d, block, item, "duplicate tag " + item.pair[0]);
      } else if (item.type == ItemType::Loop) {
        for (const LoopTag& t : item.loop.tags) {
          bool ok = names.insert(gemmi::to_lower(t.tag)).second;
          if (!ok)
            cif_fail(d, block, item, "duplicate tag " + t.tag);
        }
      } else if (item.type == ItemType::Frame) {
        bool ok = frame_names.insert(gemmi::to_lower(item.frame.name)).second;
        if (!ok)
          cif_fail(d, block, item, "duplicate save_" + item.frame.name);
      }
    }
  }
}

inline bool is_text_field(const std::string& val) {
  size_t len = val.size();
  return len > 3 && val[0] == ';' && (val[len-2] == '\n' || val[len-2] == '\r');
}

inline std::string quote(std::string v) {
  // strings with '(' happen to be quoted in mmCIF files, so we do the same
  if (v.find_first_of(" \t\r\n'\"()[]{}#") == std::string::npos &&
      v[0] != '_' && v[0] != '$' && v[0] != '\0' &&
      (v.size() > 1 || (v[0] != '.' && v[0] != '?')))
    return v;
  if (std::memchr(v.c_str(), '\n', v.size()))
    return ";" + v + "\n;";
  if (std::memchr(v.c_str(), '\'', v.size()) == nullptr)
    return "'" + v + "'";
  if (std::memchr(v.c_str(), '"', v.size()) == nullptr)
    return '"' + v + '"';
  return ";" + v + "\n;";
}

} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
