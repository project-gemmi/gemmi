// Copyright 2017 Global Phasing Ltd.
//
// struct Document that represents the CIF file (but could be also
// read from, for example, JSON file (CIF-JSON or mmJSON).

#ifndef GEMMI_CIFDOC_HPP_
#define GEMMI_CIFDOC_HPP_
#include "util.hpp"  // for starts_with
#include <algorithm> // for move, find_if, all_of, min
#include <cassert>
#include <cctype>    // for isalpha
#include <cmath>     // for NAN
#include <cstring>   // for memchr
#include <initializer_list>
#include <iosfwd>    // for size_t, istream, ptrdiff_t
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

struct TagValue {
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
  StrideIter() : cur_(nullptr), end_(nullptr), stride_(0) {}
  StrideIter(const std::vector<std::string>& vec, size_t offset,
             unsigned stride)
    : cur_(vec.data() + std::min(offset, vec.size())),
      end_(vec.data() + vec.size()),
      stride_(stride) {}
  void operator++() { cur_ = end_-cur_ > stride_ ? cur_+stride_ : end_; }
  const std::string& operator*() const { return *cur_; }
  bool operator!=(const StrideIter& other) const { return cur_ != other.cur_; }
  bool operator==(const StrideIter& other) const { return cur_ == other.cur_; }
private:
  const std::string* cur_;
  const std::string* end_;
  unsigned stride_;
};


// Loop (can by null) with position (column) corresponding to the found value.
struct LoopColumn {
  const Loop *loop;
  size_t col;
  StrideIter begin() const {
    return loop ? StrideIter(loop->values, col, loop->width())
                : StrideIter();
  }
  StrideIter end() const {
    return loop ? StrideIter(loop->values, loop->values.size(), loop->width())
                : StrideIter();
  }
};

// Some values can be given either in loop or as tag-value pairs.
// The latter case is equivalent to a loop with a single row.
// We optimized for loops, and in case of tag-values we copy the values
// into the `values` vector.
struct TableView {
  const Loop *loop;
  std::vector<int> cols;
  std::vector<std::string> values_; // used only for non-loops

  struct Row {
    const std::string* cur;
    const std::vector<int>& icols;

    const std::string& at(int n) const { return cur[icols.at(n)]; }
    const std::string& operator[](int n) const { return cur[icols[n]]; }
    size_t size() const { return icols.size(); }
    std::string str(int n) const { return as_string(at(n)); }
    struct Iter {
      const Row& parent;
      const int* cur;
      void operator++() { cur++; }
      const std::string& operator*() const { return parent.cur[*cur]; }
      bool operator!=(const Iter& other) const { return cur != other.cur; }
      bool operator==(const Iter& other) const { return cur == other.cur; }
    };
    Iter begin() const { return Iter{*this, icols.data()}; }
    Iter end() const { return Iter{*this, icols.data() + icols.size()}; }
  };

  struct Iter {
    const std::string* cur;
    const std::vector<int>* col_indices;
    size_t stride;
    void operator++() { cur += stride; }
    Row operator*() const { return Row{cur, *col_indices}; }
    bool operator!=(const Iter& other) const { return cur != other.cur; }
    bool operator==(const Iter& other) const { return cur == other.cur; }
    const std::string& get(int n) const { return cur[col_indices->at(n)]; }
  };

  bool ok() const { return loop != nullptr || !values_.empty(); }
  size_t width() const { return cols.size(); }
  size_t length() const {
    return loop ? loop->length() : (values_.empty() ? 0 : 1);
  }

  Row operator[](size_t n) const {
    const std::string* c = loop ? loop->values.data() + loop->width() * n
                                : values_.data() + values_.size() * n;
    return Row{c, cols};
  }

  Row at(size_t n) const {
    if (loop ? n >= loop->length() : n > 0 || values_.empty())
      throw std::out_of_range("No row with index " + std::to_string(n));
    return (*this)[n];
  }

  Row one() const {
    if (length() != 1)
      throw std::runtime_error("Expected one value, found " +
                                std::to_string(length()));
    return (*this)[0];
  }

  Row find_row(const std::string& s) const {
    if (loop) {
      for (size_t i = 0; i < loop->values.size(); i += loop->width())
        if (as_string(loop->values[i]) == s)
          return Row{loop->values.data() + i, cols};
    } else if (!values_.empty() && as_string(values_[0]) == s) {
      return Row{values_.data(), cols};
    }
    throw std::runtime_error("Not found in the first column: " + s);
  }

  Iter begin() const {
    if (loop)
      return Iter{loop->values.data(), &cols, loop->width()};
    if (!values_.empty())
      return Iter{values_.data(), &cols, values_.size()};
    return Iter{nullptr, nullptr, 0};
  }

  Iter end() const {
    if (loop)
      return Iter{loop->values.data() + loop->values.size(), &cols,
                  loop->width()};
    if (!values_.empty())
      return Iter{values_.data() + values_.size(), &cols, values_.size()};
    return Iter{nullptr, nullptr, 0};
  }
};

struct Item;

struct Block {
  std::string name;
  std::vector<Item> items;

  explicit Block(const std::string& name_) : name(name_) {}
  Block() {}

  // access functions
  const std::string* find_value(const std::string& tag) const;
  LoopColumn find_loop(const std::string& tag) const;
  TableView find(const std::string& prefix,
                 const std::vector<std::string>& tags) const;
  TableView find(const std::vector<std::string>& tags) const {
    return find({}, tags);
  }
  TableView find(const std::string& tag) const {
    return find({}, {tag});
  }

  // modifying functions
  void update_value(const std::string& tag, std::string v);
  bool delete_loop(const std::string& tag);
  int add_field(TableView& table, const std::string& field) const;
  Loop& clear_or_add_loop(const std::string& prefix,
                          const std::initializer_list<const char*>& tags);
  void delete_category(const std::string& prefix);
};

struct Item {
  ItemType type;
  ValueType valtype = ValueType::NotSet; // for TagValue only
  int line_number = -1;
  union {
    TagValue tv;
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
      case ItemType::Value: tv.~TagValue(); break;
      case ItemType::Loop: loop.~Loop(); break;
      case ItemType::Frame: frame.~Block(); break;
      case ItemType::Comment: tv.~TagValue(); break;
      case ItemType::Erased: break;
    }
  }

  void erase() {
    this->~Item();
    type = ItemType::Erased;
  }

  bool has_prefix(const std::string& prefix) {
    return (type == ItemType::Value && gemmi::starts_with(tv.tag, prefix)) ||
           (type == ItemType::Loop && !loop.tags.empty() &&
            gemmi::starts_with(loop.tags[0].tag, prefix));
  }

private:
  void copy_value(const Item& o) {
    if (o.type == ItemType::Value || o.type == ItemType::Comment)
      new (&tv) TagValue(o.tv);
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(o.loop);
    else if (o.type == ItemType::Frame)
      new (&frame) Block(o.frame);
  }

  void move_value(Item&& o) {
    if (o.type == ItemType::Value || o.type == ItemType::Comment)
      new (&tv) TagValue(std::move(o.tv));
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(std::move(o.loop));
    else if (o.type == ItemType::Frame)
      new (&frame) Block(std::move(o.frame));
  }
};


inline const std::string* Block::find_value(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.tv.tag == tag)
      return &i.tv.value;
  return nullptr;
}

inline void Block::update_value(const std::string& tag, std::string v) {
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

inline LoopColumn Block::find_loop(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Loop) {
      int pos = i.loop.find_tag(tag);
      if (pos != -1)
        return LoopColumn{&i.loop, static_cast<size_t>(pos)};
    }
  return LoopColumn{nullptr, 0};
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

inline TableView Block::find(const std::string& prefix,
                             const std::vector<std::string>& tags) const {
  const Loop* loop = tags.empty() ? nullptr : find_loop(prefix + tags[0]).loop;

  std::vector<int> indices;
  std::vector<std::string> values;
  if (loop) {
    indices.reserve(tags.size());
    for (const std::string& tag : tags) {
      int idx = loop->find_tag(prefix + tag);
      if (idx == -1) {
        loop = nullptr;
        indices.clear();
        break;
      }
      indices.push_back(idx);
    }
  } else {
    for (const std::string& tag : tags) {
      const std::string* v = find_value(prefix + tag);
      if (v == nullptr) {
        indices.clear();
        values.clear();
        break;
      }
      if (indices.empty()) {
        indices.reserve(tags.size());
        values.reserve(tags.size());
      }
      indices.push_back(values.size());
      values.push_back(*v);
    }
  }
  return TableView{loop, indices, values};
}

inline
int Block::add_field(TableView& table, const std::string& full_tag) const {
  if (table.loop) {
    int pos = table.loop->find_tag(full_tag);
    if (pos != -1) {
      table.cols.push_back(pos);
      return (int) table.cols.size() - 1;
    }
  } else {
    if (const std::string* v = find_value(full_tag)) {
      table.cols.push_back(table.values_.size());
      table.values_.push_back(*v);
      return (int) table.cols.size() - 1;
    }
  }
  return -1;
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
    // TODO: case-insensitive check
    bool success = names.insert(block.name).second;
    if (!success && !block.name.empty())
      throw std::runtime_error("duplicate block name: " + block.name);
  }
  // check for dups inside each block
  std::unordered_set<std::string> frame_names;
  for (const Block& block : d.blocks) {
    names.clear();
    frame_names.clear();
    for (const Item& item : block.items) {
      if (item.type == ItemType::Value) {
        bool success = names.insert(item.tv.tag).second;
        if (!success)
          cif_fail(d, block, item, "duplicate tag " + item.tv.tag);
      } else if (item.type == ItemType::Loop) {
        for (const LoopTag& t : item.loop.tags) {
          bool success = names.insert(t.tag).second;
          if (!success)
            cif_fail(d, block, item, "duplicate tag " + t.tag);
        }
      } else if (item.type == ItemType::Frame) {
        bool success = frame_names.insert(item.frame.name).second;
        if (!success)
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
