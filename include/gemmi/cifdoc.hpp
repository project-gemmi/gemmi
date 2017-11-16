// Copyright 2017 Global Phasing Ltd.
//
// struct Document that represents the CIF file (but could be also
// read from, for example, JSON file (CIF-JSON or mmJSON).

#ifndef GEMMI_CIFDOC_HPP_
#define GEMMI_CIFDOC_HPP_
#include "util.hpp"  // for starts_with, to_lower
#include <algorithm> // for move, find_if, all_of, min
#include <cassert>
#include <cstring>   // for memchr
#include <initializer_list>
#include <iosfwd>    // for size_t, ptrdiff_t
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace gemmi {
namespace cif {
using std::size_t;


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

struct Pair {
  std::string tag;
  std::string value;
};

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

  // iteration over rows
  struct Span {
    const Loop& loop;
    const std::ptrdiff_t offset;
    const std::string& at(int n) const { return loop.values.at(offset+n); }
    const std::string& operator[](int n) const { return loop.values[offset+n]; }
    size_t size() const { return loop.width(); }
    typedef std::vector<std::string>::const_iterator const_iterator;
    const_iterator begin() const { return loop.values.begin() + offset; }
    const_iterator end() const { return begin() + size(); }
  };
  struct Iter {
    const Loop& loop;
    std::vector<std::string>::const_iterator cur;
    Iter(const Loop& loop_, std::vector<std::string>::const_iterator cur_)
      : loop(loop_), cur(cur_) {}
    void operator++() { cur += loop.width(); }
    const Span operator*() const { return {loop, cur - loop.values.begin()}; }
    bool operator!=(const Iter& o) const { return cur != o.cur; }
    bool operator==(const Iter& o) const { return cur == o.cur; }
  };
  Iter begin() const { return Iter(*this, values.begin()); }
  Iter end() const { return Iter(*this, values.end()); }
};


class StrideIter {
public:
  StrideIter(const std::string* str) : cur_(str), end_(nullptr), stride_(0) {}
  StrideIter(const std::vector<std::string>& vec, size_t offset,
             unsigned stride)
    : cur_(vec.data() + std::min(offset, vec.size())),
      end_(vec.data() + vec.size()),
      stride_(stride) {}
  void operator++() {
    cur_ = (end_ && unsigned(end_-cur_) > stride_ ? cur_+stride_ : end_);
  }
  const std::string& operator*() const { return *cur_; }
  bool operator!=(const StrideIter& other) const { return cur_ != other.cur_; }
  bool operator==(const StrideIter& other) const { return cur_ == other.cur_; }
  StrideIter& to_end() { cur_ = end_; return *this; }
private:
  const std::string* cur_;
  const std::string* end_;
  unsigned stride_;
};

struct Item;
struct Block;

// Accessor to a specific loop column, or to a single value from a Pair.
struct Column {
  const Item* it;
  size_t col;  // for loop it is a column index in it->loop
  StrideIter begin() const;
  StrideIter end() const { return begin().to_end(); }
  const Loop* get_loop() const;
  const std::string* get_tag() const;
  int length() const {
    if (const Loop* loop = get_loop())
      return loop->length();
    return it ? 1 : 0;
  }
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
    struct Iter { // TODO: what should be done with absent optional tags?
      const Row& parent;
      std::vector<int>::const_iterator cur; // points into positions
      void operator++() { cur++; }
      const std::string& operator*() const { return parent.direct_at(*cur); }
      bool operator!=(const Iter& other) const { return cur != other.cur; }
      bool operator==(const Iter& other) const { return cur == other.cur; }
    };
    Iter begin() const { return Iter{*this, tab.positions.begin()}; }
    Iter end() const { return Iter{*this, tab.positions.end()}; }
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

  struct Iter {
    const Table& parent;
    size_t index;
    void operator++() { index++; }
    bool operator!=(const Iter& other) const { return index != other.index; }
    bool operator==(const Iter& other) const { return index == other.index; }
    Row operator*() const { return parent[index]; }
    const std::string& get(int n) const { return parent[index].at(n); }
  };
  Iter begin() const { return Iter{*this, 0}; }
  Iter end() const { return Iter{*this, length()}; }
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
    return pair ? &pair->value : nullptr;
  }
  Column find_loop(const std::string& tag) const;
  Column find_values(const std::string& tag) const;
  Table find(const std::string& prefix,
             const std::vector<std::string>& tags) const;
  Table find(const std::vector<std::string>& tags) const {
    return find({}, tags);
  }
  Table find(const std::string& tag) const {
    return find({}, {tag});
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
    Pair tv;
    Loop loop;
    Block frame;
  };

  explicit Item(LoopArg)
    : type{ItemType::Loop}, loop{} {}
  explicit Item(std::string&& t)
    : type{ItemType::Value}, tv{std::move(t), std::string()} {}
  Item(const std::string& t, const std::string& v)
    : type{ItemType::Value}, tv{t, v} {}
  explicit Item(FrameArg&& frame_arg)
    : type{ItemType::Frame}, frame(frame_arg.str) {}
  explicit Item(CommentArg&& comment)
    : type{ItemType::Comment}, tv{std::string(), std::move(comment.str)} {}

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
      case ItemType::Value: tv.~Pair(); break;
      case ItemType::Loop: loop.~Loop(); break;
      case ItemType::Frame: frame.~Block(); break;
      case ItemType::Comment: tv.~Pair(); break;
      case ItemType::Erased: break;
    }
  }

  void erase() {
    this->~Item();
    type = ItemType::Erased;
  }

  bool has_prefix(const std::string& prefix) const {
    return (type == ItemType::Value && gemmi::starts_with(tv.tag, prefix)) ||
           (type == ItemType::Loop && !loop.tags.empty() &&
            gemmi::starts_with(loop.tags[0].tag, prefix));
  }

private:
  void copy_value(const Item& o) {
    if (o.type == ItemType::Value || o.type == ItemType::Comment)
      new (&tv) Pair(o.tv);
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(o.loop);
    else if (o.type == ItemType::Frame)
      new (&frame) Block(o.frame);
  }

  void move_value(Item&& o) {
    if (o.type == ItemType::Value || o.type == ItemType::Comment)
      new (&tv) Pair(std::move(o.tv));
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
  return &it->tv.tag;
}

inline const Loop* Column::get_loop() const {
  return it && it->type == ItemType::Loop ? &it->loop : nullptr;
}
inline StrideIter Column::begin() const {
  if (const Loop* loop = get_loop())
    return StrideIter(loop->values, col, loop->width());
  if (it && it->type == ItemType::Value)
    return StrideIter(&it->tv.value);
  return StrideIter(nullptr);
}

inline const std::string& Table::Row::direct_at(int pos) const {
  if (pos == -1)
    throw std::out_of_range("Cannot access missing optional tag.");
  if (row_index == -1) { // tags
    if (tab.loop)
      return tab.loop->tags.at(pos).tag;
    return tab.blo.items[pos].tv.tag;
  }
  if (tab.loop)
    return tab.loop->values.at(tab.loop->width() * row_index + pos);
  return tab.blo.items[pos].tv.value;
}

inline Table::Row Table::find_row(const std::string& s) const {
  int pos = positions.at(0);
  if (loop) {
    for (size_t i = 0; i < loop->values.size(); i += loop->width())
      if (as_string(loop->values[i + pos]) == s)
        return Row{*this, static_cast<int>(i / loop->width())};
  } else if (as_string(blo.items[pos].tv.value) == s) {
    return Row{*this, 0};
  }
  throw std::runtime_error("Not found in the first column: " + s);
}

inline const Item* Block::find_pair_item(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.tv.tag == tag)
      return &i;
  return nullptr;
}

inline const Pair* Block::find_pair(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.tv.tag == tag)
      return &i.tv;
  return nullptr;
}

inline void Block::set_pair(const std::string& tag, std::string v) {
  if (tag[0] != '_')
    throw std::runtime_error("Tag should start with '_', got: " + tag);
  for (Item& i : items) {
    if (i.type == ItemType::Value && i.tv.tag == tag) {
      i.tv.value = v;
      return;
    }
    if (i.type == ItemType::Loop && i.loop.find_tag(tag) != -1) {
      i.loop.~Loop();
      i.type = ItemType::Value;
      i.tv = { tag, v };
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
      if (i.tv.tag == tag)
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
      tag = &item.tv.tag;
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
        bool ok = names.insert(gemmi::to_lower(item.tv.tag)).second;
        if (!ok)
          cif_fail(d, block, item, "duplicate tag " + item.tv.tag);
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
