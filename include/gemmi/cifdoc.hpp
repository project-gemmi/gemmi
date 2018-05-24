// Copyright 2017 Global Phasing Ltd.
//
// struct Document that represents the CIF file (but could be also
// read from, for example, JSON file (CIF-JSON or mmJSON).

#ifndef GEMMI_CIFDOC_HPP_
#define GEMMI_CIFDOC_HPP_
#include "util.hpp"  // for starts_with, to_lower, fail
#include "iterator.hpp"  // for StrideIter, IndirectIter
#include <algorithm> // for move, find_if, all_of, min
#include <array>
#include <cstring>   // for memchr
#include <initializer_list>
#include <iosfwd>    // for size_t, ptrdiff_t
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#if defined(_MSC_VER)
#pragma warning(push)
// warning C4244: an integer type is converted to a smaller integer type
#pragma warning(disable: 4244)
// warning C4267: conversion from 'size_t' to 'type', possible loss of data
#pragma warning(disable: 4267)
#endif

namespace gemmi {
namespace cif {
using std::size_t;
using gemmi::fail;

enum class ItemType : unsigned char {
  Pair,
  Loop,
  Frame,
  Comment,
  Erased,
};

inline void assert_tag(const std::string& tag) {
  if (tag[0] != '_')
    fail("Tag should start with '_', got: " + tag);
}

inline void ensure_mmcif_category(std::string& cat) {
  if (cat[0] != '_')
    fail("Category should start with '_', got: " + cat);
  if (*(cat.end() - 1) != '.')
    cat += '.';
}

inline bool is_null(const std::string& value) {
  return value.size() == 1 && (value[0] == '?' || value[0] == '.');
}

inline std::string as_string(const std::string& value) {
  if (value.empty() || is_null(value))
    return "";
  if (value[0] == '"' || value[0] == '\'')
    return std::string(value.begin() + 1, value.end() - 1);
  if (value[0] == ';' && value.size() > 2 && *(value.end() - 2) == '\n') {
    bool crlf = *(value.end() - 3) == '\r';
    return std::string(value.begin() + 1, value.end() - (crlf ? 3 : 2));
  }
  return value;
}

inline std::string as_string(const std::string* value) {
  return value ? as_string(*value) : std::string();
}

inline char as_char(const std::string& value, char null) {
  if (is_null(value))
    return null;
  if (value.size() < 2)
    return value[0];
  const std::string s = as_string(value);
  if (s.size() < 2)
    return s[0];
  fail("Not a single character: " + value);
}

using Pair = std::array<std::string, 2>;

// used only as arguments when creating Item
struct LoopArg {};
struct FrameArg { std::string str; };
struct CommentArg { std::string str; };

struct Loop {
  std::vector<std::string> tags;
  std::vector<std::string> values;

  // search and access
  int find_tag(const std::string& tag) const {
    auto f = std::find_if(tags.begin(), tags.end(),
                          [&tag](const std::string& t) { return t == tag; });
    return f == tags.end() ? -1 : f - tags.begin();
  }
  size_t width() const { return tags.size(); }
  size_t length() const { return values.size() / tags.size(); }
  const std::string& val(size_t row, size_t col) const {
    return values[row * tags.size() + col];
  }

  void clear() { tags.clear(); values.clear(); }

  template <typename T> void add_row(T new_values, int pos=-1) {
    if (new_values.size() != tags.size())
      fail("add_row(): wrong row length.");
    auto it = values.end();
    if (pos >= 0 && pos * width() < values.size())
      it = values.begin() + pos * tags.size();
    values.insert(it, new_values.begin(), new_values.end());
  }
  void add_row(std::initializer_list<std::string> new_values, int pos=-1) {
    add_row<std::initializer_list<std::string>>(new_values, pos);
  }

  void set_all_values(std::vector<std::vector<std::string>> columns);
};


struct Item;
struct Block;

// Accessor to a specific loop column, or to a single value from a Pair.
class Column {
public:
  Column() : item_(nullptr) {}
  Column(Item* item, size_t col) : item_(item), col_(col) {}
  using iterator = StrideIter<std::string>;
  iterator begin();
  iterator end();
  using const_iterator = StrideIter<const std::string>;
  const_iterator begin() const { return const_cast<Column*>(this)->begin(); }
  const_iterator end() const { return const_cast<Column*>(this)->end(); }

  Loop* get_loop() const;
  std::string* get_tag();
  const std::string* get_tag() const {
    return const_cast<Column*>(this)->get_tag();
  }
  int length() const {
    if (const Loop* loop = get_loop())
      return loop->length();
    return item_ ? 1 : 0;
  }
  explicit operator bool() const { return item_ != nullptr ; }
  std::string& operator[](int n);
  std::string& at(int n) {
    if (n < 0)
      n += length();
    if (n < 0 || n >= length())
      throw std::out_of_range("Cannot access element " + std::to_string(n) +
          " in Column with length " + std::to_string(length()));
    return operator[](n);
  }
  const std::string& at(int n) const {
    return const_cast<Column*>(this)->at(n);
  }
  std::string str(int n) const { return as_string(at(n)); }
  const Item* item() const { return item_; }
  Item* item() { return item_; }
  size_t col() const { return col_; }

private:
  Item* item_;
  size_t col_;  // for loop this is a column index in item_->loop
};

// Some values can be given either in loop or as tag-value pairs.
// The latter case is equivalent to a loop with a single row.
// We optimized for loops, and in case of tag-values we copy the values
// into the `values` vector.
struct Table {
  Item* loop_item;
  Block& bloc;
  std::vector<int> positions;

  struct Row {
    Table& tab;
    int row_index;

    std::string& value_at(int pos);
    const std::string& value_at(int pos) const {
      return const_cast<Row*>(this)->value_at(pos);
    }
    std::string& at(int n) {
      return value_at(tab.positions.at(n < 0 ? n + size() : n));
    }
    const std::string& at(int n) const { return const_cast<Row*>(this)->at(n); }
    std::string& operator[](int n);
    const std::string& operator[](int n) const {
      return const_cast<Row*>(this)->operator[](n);
    }
    std::string* ptr_at(int n) {
      int pos = tab.positions.at(n < 0 ? n + size() : n);
      return pos >= 0 ? &value_at(pos) : nullptr;
    }
    const std::string* ptr_at(int n) const {
      return const_cast<Row*>(this)->ptr_at(n);
    }
    bool has(int n) const { return tab.has_column(n); }
    bool has2(int n) const { return has(n) && !cif::is_null(operator[](n)); }
    size_t size() const { return tab.width(); }
    std::string str(int n) const { return as_string(at(n)); }
    using iterator = IndirectIter<Row, std::string>;
    using const_iterator = IndirectIter<const Row, const std::string>;
    iterator begin() { return iterator({this, tab.positions.begin()}); }
    iterator end() { return iterator({this, tab.positions.end()}); }
    const_iterator begin() const {
      return const_iterator({this, tab.positions.begin()});
    }
    const_iterator end() const {
      return const_iterator({this, tab.positions.end()});
    }
  };

  bool ok() const { return !positions.empty(); }
  size_t width() const { return positions.size(); }
  size_t length() const;
  bool has_column(int n) const { return positions.at(n) >= 0; }
  Row tags() { return Row{*this, -1}; }
  Row operator[](int n) { return Row{*this, n}; }

  Row at(int n) {
    if (n < 0)
      n += length();
    if (n < 0 || static_cast<size_t>(n) >= length())
      throw std::out_of_range("No row with index " + std::to_string(n));
    return (*this)[n];
  }

  Row one() {
    if (length() != 1)
      fail("Expected one value, found " + std::to_string(length()));
    return (*this)[0];
  }

  Row find_row(const std::string& s);

  Column column(int n);

  // tries exact tag match first, looks for suffix later
  Column find_column(const std::string& suffix) {
    int w = width();
    Row tag_row = tags();
    for (int i = 0; i != w; ++i)
      if (tag_row[i] == suffix)
        return column(i);
    for (int i = 0; i != w; ++i)
      if (gemmi::ends_with(tag_row[i], suffix))
        return column(i);
    fail("Column name or suffix not found: " + suffix);
  }

  void erase();

  // It is not a proper input iterator, but just enough for using range-for.
  struct iterator {
    Table& parent;
    size_t index;
    void operator++() { index++; }
    bool operator==(const iterator& o) const { return index == o.index; }
    bool operator!=(const iterator& o) const { return index != o.index; }
    Row operator*() { return parent[index]; }
    const std::string& get(int n) const { return parent[index].at(n); }
  };
  iterator begin() { return iterator{*this, 0}; }
  iterator end() { return iterator{*this, length()}; }
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
  Column find_loop(const std::string& tag);
  Column find_values(const std::string& tag);
  Table find(const std::string& prefix,
             const std::vector<std::string>& tags);
  Table find(const std::vector<std::string>& tags) { return find({}, tags); }
  Table find_any(const std::string& prefix,
                 const std::vector<std::string>& tags);

  // modifying functions
  void set_pair(const std::string& tag, std::string v);

  Loop& init_loop(const std::string& prefix, std::vector<std::string> tags) {
    return setup_loop(find_any(prefix, tags), prefix, std::move(tags));
  }

  // mmCIF specific functions
  std::vector<std::string> get_mmcif_category_names() const;
  Table find_mmcif_category(std::string cat);

  Loop& init_mmcif_loop(std::string cat, std::vector<std::string> tags) {
    ensure_mmcif_category(cat);
    return setup_loop(find_mmcif_category(cat), cat, std::move(tags));
  }

private:
  Loop& setup_loop(Table&& tab, const std::string& prefix,
                   std::vector<std::string>&& tags);
};

struct Item {
  ItemType type;
  int line_number = -1;
  union {
    Pair pair;
    Loop loop;
    Block frame;
  };

  explicit Item(LoopArg)
    : type{ItemType::Loop}, loop{} {}
  explicit Item(std::string&& t)
    : type{ItemType::Pair}, pair{{std::move(t), std::string()}} {}
  Item(const std::string& t, const std::string& v)
    : type{ItemType::Pair}, pair{{t, v}} {}
  explicit Item(FrameArg&& frame_arg)
    : type{ItemType::Frame}, frame(frame_arg.str) {}
  explicit Item(CommentArg&& comment)
    : type{ItemType::Comment}, pair{{std::string(), std::move(comment.str)}} {}

  Item(Item&& o) noexcept
      : type(o.type), line_number(o.line_number) {
    move_value(std::move(o));
  }
  Item(const Item&& o)
      : type(o.type), line_number(o.line_number) {
    copy_value(o);
  }
  Item(const Item& o)
      : type(o.type), line_number(o.line_number) {
    copy_value(o);
  }
  Item(Item& o) : Item(static_cast<const Item&>(o)) {}

  ~Item() {
    switch (type) {
      case ItemType::Pair: pair.~Pair(); break;
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
    return (type == ItemType::Pair && gemmi::starts_with(pair[0], prefix)) ||
           (type == ItemType::Loop && !loop.tags.empty() &&
            gemmi::starts_with(loop.tags[0], prefix));
  }

  void set_value(Item&& o) {
    if (type == o.type) {
      switch (type) {
        case ItemType::Pair: pair = std::move(o.pair); break;
        case ItemType::Loop: loop = std::move(o.loop); break;
        case ItemType::Frame: frame = std::move(o.frame); break;
        case ItemType::Comment: pair = std::move(o.pair); break;
        case ItemType::Erased: break;
      }
    } else {
      this->~Item();
      type = o.type;
      move_value(std::move(o));
    }
  }

private:
  void copy_value(const Item& o) {
    if (o.type == ItemType::Pair || o.type == ItemType::Comment)
      new (&pair) Pair(o.pair);
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(o.loop);
    else if (o.type == ItemType::Frame)
      new (&frame) Block(o.frame);
  }

  void move_value(Item&& o) {
    if (o.type == ItemType::Pair || o.type == ItemType::Comment)
      new (&pair) Pair(std::move(o.pair));
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(std::move(o.loop));
    else if (o.type == ItemType::Frame)
      new (&frame) Block(std::move(o.frame));
  }
};

inline void Loop::set_all_values(std::vector<std::vector<std::string>> columns){
  size_t w = columns.size();
  if (w != width())
    fail("set_all_values(): expected " + std::to_string(width()) +
         " columns, got " + std::to_string(w));
  if (w == 0)
    return;
  size_t h = columns[0].size();
  for (auto& col : columns)
    if (col.size() != h)
      fail("set_all_values(): all columns must have the same length");
  values.resize(w * h);
  for (size_t i = 0; i != h; ++i)
    for (size_t j = 0; j != w; ++j)
      values[w * i + j] = std::move(columns[j][i]);
}

inline std::string* Column::get_tag() {
  if (!item_)
    return nullptr;
  if (Loop* loop = get_loop())
    return &loop->tags.at(col_);
  return &item_->pair[0];
}

inline Loop* Column::get_loop() const {
  return item_ && item_->type == ItemType::Loop ? &item_->loop : nullptr;
}
inline Column::iterator Column::begin() {
  if (Loop* loop = get_loop())
    return iterator({loop->values.data(), col_, (unsigned) loop->width()});
  if (item_ && item_->type == ItemType::Pair)
    return iterator({&item_->pair[1], 0, 1});
  return iterator();
}

inline Column::iterator Column::end() {
  if (Loop* loop = get_loop())
    return iterator({loop->values.data() + loop->values.size(),
                    col_, (unsigned) loop->width()});
  if (item_ && item_->type == ItemType::Pair)
    return iterator({&item_->pair[1] + 1, 0, 1});
  return iterator();
}

inline std::string& Column::operator[](int n) {
  if (Loop* loop = get_loop())
    return loop->values[n * loop->width() + col_];
  return item_->pair[1];
}

inline std::string& Table::Row::operator[](int n) {
  int pos = tab.positions[n];
  if (tab.loop_item) {
    Loop& loop = tab.loop_item->loop;
    if (row_index == -1) // tags
      return loop.tags[pos];
    return loop.values[loop.width() * row_index + pos];
  }
  return tab.bloc.items[pos].pair[row_index == -1 ? 0 : 1];
}

inline std::string& Table::Row::value_at(int pos) {
  if (pos == -1)
    throw std::out_of_range("Cannot access missing optional tag.");
  if (row_index == -1) { // tags
    if (tab.loop_item)
      return tab.loop_item->loop.tags.at(pos);
    return tab.bloc.items[pos].pair[0];
  }
  if (tab.loop_item) {
    Loop& loop = tab.loop_item->loop;
    return loop.values.at(loop.width() * row_index + pos);
  }
  return tab.bloc.items[pos].pair[1];
}

inline size_t Table::length() const {
  return loop_item ? loop_item->loop.length() : (positions.empty() ? 0 : 1);
}

inline Table::Row Table::find_row(const std::string& s) {
  int pos = positions.at(0);
  if (loop_item) {
    const Loop& loop = loop_item->loop;
    for (size_t i = 0; i < loop.values.size(); i += loop.width())
      if (as_string(loop.values[i + pos]) == s)
        return Row{*this, static_cast<int>(i / loop.width())};
  } else if (as_string(bloc.items[pos].pair[1]) == s) {
    return Row{*this, 0};
  }
  fail("Not found in the first column: " + s);
}

inline Column Table::column(int n) {
  int pos = positions.at(n);
  if (pos == -1)
    fail("Cannot access absent column");
  if (loop_item)
    return Column(loop_item, pos);
  return Column(&bloc.items[pos], 0);
}

inline void Table::erase() {
  if (loop_item)
    loop_item->erase();
  else
    for (int pos : positions)
      bloc.items[pos].erase();
}

inline const Item* Block::find_pair_item(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Pair && i.pair[0] == tag)
      return &i;
  return nullptr;
}

inline const Pair* Block::find_pair(const std::string& tag) const {
  const Item* item = find_pair_item(tag);
  return item ? &item->pair : nullptr;
}

inline void Block::set_pair(const std::string& tag, std::string v) {
  assert_tag(tag);
  for (Item& i : items) {
    if (i.type == ItemType::Pair && i.pair[0] == tag) {
      i.pair[1] = v;
      return;
    }
    if (i.type == ItemType::Loop && i.loop.find_tag(tag) != -1) {
      i.set_value(Item(tag, v));
      return;
    }
  }
  items.emplace_back(tag, v);
}

inline Column Block::find_loop(const std::string& tag) {
  Column c = find_values(tag);
  return c.item() && c.item()->type == ItemType::Loop ? c : Column();
}

inline Column Block::find_values(const std::string& tag) {
  for (Item& i : items)
    if (i.type == ItemType::Loop) {
      int pos = i.loop.find_tag(tag);
      if (pos != -1)
        return Column{&i, static_cast<size_t>(pos)};
    } else if (i.type == ItemType::Pair) {
      if (i.pair[0] == tag)
        return Column{&i, 0};
    }
  return Column{nullptr, 0};
}

inline std::vector<std::string> Block::get_mmcif_category_names() const {
  std::vector<std::string> cats;
  for (const Item& item : items) {
    const std::string* tag = nullptr;
    if (item.type == ItemType::Pair)
      tag = &item.pair[0];
    else if (item.type == ItemType::Loop && !item.loop.tags.empty())
      tag = &item.loop.tags[0];
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

inline Loop& Block::setup_loop(Table&& tab, const std::string& prefix,
                               std::vector<std::string>&& tags) {
  Item *item;
  if (tab.loop_item) {
    item = tab.loop_item;
    item->loop.clear();
  } else if (tab.ok()) {
    item = &tab.bloc.items.at(tab.positions[0]);
    tab.erase();
    item->set_value(Item(LoopArg{}));
  } else {
    items.emplace_back(LoopArg{});
    item = &items.back();
  }
  for (std::string& tag : tags) {
    tag.insert(0, prefix);
    assert_tag(tag);
  }
  item->loop.tags = std::move(tags);
  return item->loop;
}

inline Table Block::find(const std::string& prefix,
                         const std::vector<std::string>& tags) {
  Item* loop_item = nullptr;
  if (!tags.empty())
    loop_item = find_loop(prefix + tags[0]).item();

  std::vector<int> indices;
  indices.reserve(tags.size());
  if (loop_item) {
    for (const std::string& tag : tags) {
      std::string full_tag = prefix + (tag[0] != '?' ? tag : tag.substr(1));
      int idx = loop_item->loop.find_tag(full_tag);
      if (idx == -1 && tag[0] != '?') {
        loop_item = nullptr;
        indices.clear();
        break;
      }
      indices.push_back(idx);
    }
  } else {
    for (const std::string& tag : tags) {
      std::string full_tag = prefix + (tag[0] != '?' ? tag : tag.substr(1));
      if (const Item* p = find_pair_item(full_tag)) {
        indices.push_back(p - items.data());
      } else if (tag[0] == '?') {
        indices.push_back(-1);
      } else {
        indices.clear();
        break;
      }
    }
  }
  return Table{loop_item, *this, indices};
}

inline Table Block::find_any(const std::string& prefix,
                             const std::vector<std::string>& tags) {
  std::vector<int> indices;
  for (auto tag = tags.begin(); tag != tags.end(); ++tag) {
    Column column = find_values(prefix + *tag);
    if (Item* item = column.item()) {
      if (item->type == ItemType::Loop) {
        indices.push_back(column.col());
        while (++tag != tags.end()) {
          int idx = item->loop.find_tag(prefix + *tag);
          if (idx != -1)
            indices.push_back(idx);
        }
        return Table{item, *this, indices};
      } else {
        indices.push_back(item - items.data());
        while (++tag != tags.end())
          if (const Item* p = find_pair_item(prefix + *tag))
            indices.push_back(p - items.data());
        return Table{nullptr, *this, indices};
      }
    }
  }
  return Table{nullptr, *this, indices};
}

inline Table Block::find_mmcif_category(std::string cat) {
  ensure_mmcif_category(cat);
  std::vector<int> indices;
  for (Item& i : items)
    if (i.has_prefix(cat)) {
      if (i.type == ItemType::Loop) {
        indices.clear();
        indices.resize(i.loop.tags.size());
        for (size_t j = 0; j != indices.size(); ++j) {
          indices[j] = j;
          const std::string& tag = i.loop.tags[j];
          if (!starts_with(tag, cat))
            fail("Tag " + tag + " in loop with " + cat);
        }
        return Table{&i, *this, indices};
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

  Block& add_new_block(const std::string& name, int pos=-1) {
    if (find_block(name))
      fail("Block with such name already exists: " + name);
    if (pos > 0 && static_cast<size_t>(pos) > blocks.size())
      throw std::out_of_range("add_new_block(): invalid position");
    return *blocks.emplace(pos < 0 ? blocks.end() : blocks.begin() + pos, name);
  }

  void clear() noexcept {
    source.clear();
    blocks.clear();
    items_ = nullptr;
  }

  // returns blocks[0] if the document has exactly one block (like mmCIF)
  Block& sole_block() {
    if (blocks.size() > 1)
      fail("single data block expected, got " + std::to_string(blocks.size()));
    return blocks.at(0);
  }
  const Block& sole_block() const {
    return const_cast<Document*>(this)->sole_block();
  }

  Block* find_block(const std::string& name) {
    for (Block& b : blocks)
      if (b.name == name)
        return &b;
    return nullptr;
  }
  const Block* find_block(const std::string& name) const {
    return const_cast<Document*>(this)->find_block(name);
  }
};


[[noreturn]]
inline void cif_fail(const Document& d, const Block& b, const Item& item,
                     const std::string& s) {
  fail(d.source + ":" + std::to_string(item.line_number) +
       " in data_" + b.name + ": " + s);
}

// Throw an error if any block name, frame name or tag is duplicated.
inline void check_duplicates(const Document& d) {
  // check for duplicate block names (except empty "" which is global_)
  std::unordered_set<std::string> names;
  for (const Block& block : d.blocks) {
    bool ok = names.insert(gemmi::to_lower(block.name)).second;
    if (!ok && !block.name.empty())
      fail("duplicate block name: " + block.name);
  }
  // check for dups inside each block
  std::unordered_set<std::string> frame_names;
  for (const Block& block : d.blocks) {
    names.clear();
    frame_names.clear();
    for (const Item& item : block.items) {
      if (item.type == ItemType::Pair) {
        bool ok = names.insert(gemmi::to_lower(item.pair[0])).second;
        if (!ok)
          cif_fail(d, block, item, "duplicate tag " + item.pair[0]);
      } else if (item.type == ItemType::Loop) {
        for (const std::string& t : item.loop.tags) {
          bool ok = names.insert(gemmi::to_lower(t)).second;
          if (!ok)
            cif_fail(d, block, item, "duplicate tag " + t);
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

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
