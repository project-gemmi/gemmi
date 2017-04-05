// Copyright 2017 Global Phasing Ltd.
//
// CIF parser (based on PEGTL), struct Document that represents the CIF file,
// and a set of actions for the parser to prepare Document.

#ifndef GEMMI_CIF_HH_
#define GEMMI_CIF_HH_
#include "numb.hh" // is_int() for infer_valtypes()
#include <cassert>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <new>
#include <unordered_set>

#include <tao/pegtl.hpp>
#ifdef CIF_VALIDATE_SHOW_TRACE
#include <tao/pegtl/trace.hpp>
#endif

namespace gemmi {
namespace cif {
using std::size_t;
namespace pegtl = tao::pegtl;


// **** grammar rules, named similarly as in the CIF 1.1 spec ****
namespace rules {

  using namespace pegtl;

  // (letter) refers to sections in Table 2.2.7.1 in Vol.G of ITfC (2006).

  // (g) Character sets.
  // OrdinaryCharacter: ! % &  ()*+,-./0-9:  <=>?@A-Z[]  \ ^  `a-z{|}~

  // !"#$%&'()*+,-./0-9:;<=>?@A-Z[\]^_`a-z{|}~
  struct nonblank_ch : range<'!', '~'> {};

  // ascii space is just before '!'
  struct anyprint_ch : ranges<' ', '~', '\t'> {};


  // (f) White space and comments.
  struct comment : if_must<one<'#'>, until<eolf>> {};
  struct whitespace : plus<sor<one<' ','\n','\r','\t'>, comment>> {};

  // (b) Reserved words.
  struct str_data : TAOCPP_PEGTL_ISTRING("data_") {};
  struct str_loop : TAOCPP_PEGTL_ISTRING("loop_") {};
  struct str_global : TAOCPP_PEGTL_ISTRING("global_") {};
  struct str_save : TAOCPP_PEGTL_ISTRING("save_") {};
  struct str_stop : TAOCPP_PEGTL_ISTRING("stop_") {};
  struct keyword : sor<str_data, str_loop, str_global, str_save, str_stop> {};

  // (e) Character strings and text fields.
  template <typename Q>
  struct endq : seq<Q, at<sor<one<' ','\n','\r','\t','#'>, eof>>> {};
  template <typename Q> struct quoted_tail : until<endq<Q>, anyprint_ch> {};
  template <typename Q> struct quoted : if_must<Q, quoted_tail<Q>> {};
  struct singlequoted : quoted<one<'\''>> {};
  struct doublequoted : quoted<one<'"'>> {};
  struct field_sep : seq<bol, one<';'>> {};
  // CIF 2.0 requires whitespace after text field, so it'd be:
  // until<endq<field_sep>> instead of until<field_sep>.
  struct textfield : if_must<field_sep, until<field_sep>> {};
  struct unquoted : seq<not_at<keyword>, not_at<one<'_','$','#'>>,
                        plus<nonblank_ch>> {};

  // (c) Tags and values.

  // (a) Basic structure of CIF.
  struct datablockname : plus<nonblank_ch> {};
  struct datablockheading : sor<if_must<str_data, datablockname>,
                                str_global> {};
  struct tag : seq<one<'_'>, plus<nonblank_ch>> {};
  // simple unquoted value - for a typical mmCIF file it is faster
  // to check first this risking backtracking.
  struct simunq : seq<plus<ranges<'(', ':', '<', 'Z', 'a', 'z'>>,
                      at<one<' ','\n','\r','\t'>>> {};
  struct value: sor<simunq, singlequoted, doublequoted, textfield, unquoted> {};
  struct loop_tag : tag {};
  struct loop_value : value {};
  struct loop: if_must<str_loop, plus<seq<whitespace, loop_tag, discard>>,
                                 star<seq<whitespace, loop_value, discard>>,
                                 opt<seq<whitespace, str_stop>>> {};
  struct dataitem: if_must<tag, whitespace, value, discard> {};
  struct framename : plus<nonblank_ch> {};
  struct endframe : str_save {};
  struct frame : if_must<str_save, framename, whitespace,
                         star<sor<dataitem, loop>, whitespace>,
                         endframe> {};
  struct datablock : seq<datablockheading,
                         star<whitespace, sor<dataitem, loop, frame>>> {};
  struct file : must<opt<whitespace>,
                     until<eof, list_tail<datablock, whitespace>>> {};

} // namespace rules


// **** error messages ****

template<typename Rule> struct Errors : public pegtl::normal<Rule> {
  static const std::string msg;

  template<typename Input, typename ... States>
  static void raise(const Input& in, States&& ...) {
    throw pegtl::parse_error(msg, in);
  }
};
#define error_msg(x) template<> const std::string Errors<x>::msg
// TODO: error "expected data_ keyword
error_msg(rules::quoted_tail<rules::one<'\''>>) = "unterminated 'string'";
error_msg(rules::quoted_tail<rules::one<'"'>>) = "unterminated \"string\"";
error_msg(rules::quoted_tail<rules::field_sep>) = "unterminated text field";
error_msg(rules::value) = "expected value";
error_msg(rules::datablockname) = "unnamed DataBlock";
error_msg(rules::framename) = "unnamed save_ frame";
#undef error_msg
template<typename T> const std::string Errors<T>::msg = "parse error";


// **** data storage that preserves the order, comments, etc ****

enum class ItemType : unsigned char {
  Value,
  Loop,
  Frame,
};

enum class ValueType : unsigned char {
  NotSet,
  Char, // Line/Text?
  Numb, // Int/Float?
  Dot,
  QuestionMark,
};

std::string value_type_to_str(ValueType v) {
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
  if (value.empty())
    return value;
  else if (value[0] == '"' || value[0] == '\'')
    return std::string(value.begin() + 1, value.end() - 1);
  else if (value[0] == ';' && value.size() > 2 &&
           *(value.end() - 2) == '\n')
    return std::string(value.begin() + 1, value.end() - 2);
  else
    return value;
}

struct TagValue {
  std::string tag;
  std::string value;
};

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
};


class StrideIter {
public:
  StrideIter() : cur_(nullptr), end_(nullptr), stride_(0) {}
  StrideIter(const std::vector<std::string>& vec, size_t offset, int stride)
    : cur_(vec.data() + std::min(offset, vec.size())),
      end_(vec.data() + vec.size()),
      stride_(stride) {}
  void operator++() { cur_ = end_-cur_ > stride_ ? cur_+stride_ : end_; }
  const std::string& operator*() const { return *cur_; };
  bool operator!=(const StrideIter& other) const { return cur_ != other.cur_; }
  bool operator==(const StrideIter& other) const { return cur_ == other.cur_; }
private:
  const std::string* cur_;
  const std::string* end_;
  unsigned int stride_;
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

struct LoopTable {
  const Loop *loop;
  std::vector<int> cols;
  struct Row {
    const std::string* cur;
    const std::vector<int>& icols;
    const std::string& raw(int n) const { return cur[icols.at(n)]; }
    const std::string& operator[](int n) const { return raw(n); }
    std::string as_str(int n) const { return as_string(raw(n)); }
    double as_num(int n) const { return as_number(raw(n)); }
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
  };
  Iter begin() const {
    return loop ? Iter{loop->values.data(), &cols, loop->width()}
                : Iter{nullptr, nullptr, 0};
  }
  Iter end() const {
    return loop ? Iter{loop->values.data() + loop->values.size(), &cols,
                       loop->width()}
                : Iter{nullptr, nullptr, 0};
  }
};

struct Item;

struct Block {
  std::string name;
  std::vector<Item> items;

  explicit Block(const std::string& name_) : name(name_) {}
  Block() {}

  const std::string* find_value(const std::string& tag) const;
  LoopColumn find_loop(const std::string& tag) const;
  LoopTable find_loop_values(const std::string& prefix,
                             const std::vector<std::string>& tags) const;
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

  explicit Item(int)
    : type{ItemType::Loop}, loop{} {}
  explicit Item(std::string&& t)
    : type{ItemType::Value}, tv{std::move(t), std::string()} {}
  Item(const std::string& t, const std::string& v)
    : type{ItemType::Value}, tv{t, v} {}
  Item(const std::string& s, int)
    : type{ItemType::Frame}, frame(s) {}

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
  Item(Item& o) : Item(static_cast<const Item &>(o)) {}

  ~Item() {
    switch (type) {
      case ItemType::Value: tv.~TagValue(); break;
      case ItemType::Loop: loop.~Loop(); break;
      case ItemType::Frame: frame.~Block(); break;
    }
  }

private:
  void copy_value(const Item& o) {
    if (o.type == ItemType::Value)
      new (&tv) TagValue(o.tv);
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(o.loop);
    else if (o.type == ItemType::Frame)
      new (&frame) Block(o.frame);
  }

  void move_value(Item&& o) {
    if (o.type == ItemType::Value)
      new (&tv) TagValue(std::move(o.tv));
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(std::move(o.loop));
    else if (o.type == ItemType::Frame)
      new (&frame) Block(std::move(o.frame));
  }
};


struct Comment {
  int line_number;
  std::string text;
  Comment(int line, std::string s) : line_number(line), text(s) {}
};


inline const std::string* Block::find_value(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.tv.tag == tag)
      return &i.tv.value;
  return nullptr;
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

inline LoopTable Block::find_loop_values(const std::string& prefix,
    const std::vector<std::string>& tags) const {
  if (tags.empty())
    return LoopTable{nullptr, {}};
  LoopColumn c0 = find_loop(prefix + tags[0]);
  if (!c0.loop)
    return LoopTable{nullptr, {}};
  std::vector<int> indices;
  indices.reserve(tags.size());
  for (const std::string& tag : tags) {
    int idx = c0.loop->find_tag(prefix + tag);
    if (idx == -1)
      return LoopTable{nullptr, {}};
    indices.push_back(idx);
  }
  return LoopTable{c0.loop, indices};
}

struct Document {
  Document() : items_{nullptr} {}
  explicit Document(const std::string& path) : items_{nullptr} {
    read_file(path);
  }

  void read_file(const std::string& filename);
  void read_string(const std::string& data, const std::string& name="string");
  void read_memory(const char* data, const size_t size, const char* name);
  void read_cstream(std::FILE *f, const char* name, size_t maximum);
  void read_istream(std::istream &is, const char* name, size_t maximum);

  // returns blocks[0] if the document has exactly one block (like mmCIF)
  const Block& sole_block() const {
    if (blocks.size() > 1)
      throw std::runtime_error(std::to_string(blocks.size()) + " data blocks,"
                               " a single block was expected");
    return blocks.at(0);
  }

  std::string source;
  std::vector<Block> blocks;
  std::vector<Comment> comments;

  // implementation detail
  std::vector<Item>* items_; // items of the currently parsed block or frame

private:
  void after_read(const std::string& name);
};

// **** parsing actions that fill the storage ****

template<typename Rule> struct Action : pegtl::nothing<Rule> {};

template<> struct Action<rules::datablockname> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.blocks.emplace_back(in.string());
    out.items_ = &out.blocks.back().items;
  }
};
template<> struct Action<rules::str_global> {
  template<typename Input> static void apply(const Input&, Document& out) {
    out.blocks.emplace_back();
    out.items_ = &out.blocks.back().items;
  }
};
template<> struct Action<rules::framename> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.items_->emplace_back(in.string(), 0);
    out.items_->back().line_number = in.line();
    out.items_ = &out.items_->back().frame.items;
  }
};
template<> struct Action<rules::endframe> {
  template<typename Input> static void apply(const Input&, Document& out) {
    out.items_ = &out.blocks.back().items;
  }
};
template<> struct Action<rules::tag> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.items_->emplace_back(in.string());
    out.items_->back().line_number = in.line();
  }
};
template<> struct Action<rules::value> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items_->back();
    assert(last_item.type == ItemType::Value);
    last_item.tv.value = in.string();
  }
};
template<> struct Action<rules::str_loop> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.items_->emplace_back(0);
    out.items_->back().line_number = in.line();
  }
};
template<> struct Action<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items_->back();
    assert(last_item.type == ItemType::Loop);
    last_item.loop.tags.emplace_back(in.string());
    last_item.loop.tags.back().line_number = in.line();
  }
};
template<> struct Action<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items_->back();
    assert(last_item.type == ItemType::Loop);
    last_item.loop.values.emplace_back(in.string());
  }
};
template<> struct Action<rules::loop> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items_->back();
    assert(last_item.type == ItemType::Loop);
    const Loop& loop = last_item.loop;
    if (loop.values.size() % loop.tags.size() != 0)
      throw pegtl::parse_error("Wrong number of values in the loop", in);
  }
};
template<> struct Action<rules::comment> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    // FIXME: should we ignore silly empty non-comments in mmCIFs from PDB?
    out.comments.emplace_back(in.line(), in.string());
  }
};


void throw_validation_err(const Document& d, const Block& b, const Item& item,
                          const std::string& s) {
  throw std::runtime_error(d.source + ":" + std::to_string(item.line_number) +
                           " in data_" + b.name + ": " + s);
}

// Throw an error if any block name, frame name or tag is duplicated.
inline void check_duplicates(const Document& d) {
  // check for duplicate block names (except empty "" which is global_)
  std::unordered_set<std::string> names;
  for (const Block& block : d.blocks) {
    if (block.name.empty())
      continue;
    bool success = names.insert(block.name).second;
    if (!success)
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
          throw_validation_err(d, block, item, "duplicate tag " + item.tv.tag);
      } else if (item.type == ItemType::Loop) {
        for (const LoopTag& t : item.loop.tags) {
          bool success = names.insert(t.tag).second;
          if (!success)
            throw_validation_err(d, block, item, "duplicate tag " + t.tag);
        }
      } else if (item.type == ItemType::Frame) {
        bool success = frame_names.insert(item.frame.name).second;
        if (!success)
          throw_validation_err(d, block, item,
                               "duplicate save_" + item.frame.name);
      }
    }
  }
}

inline void Document::read_file(const std::string& filename) {
  pegtl::parse_file<rules::file, Action, Errors>(filename, *this);
  //pegtl::read_parser(filename).parse<rules::file, Action, Errors>(*this);
  after_read(filename);
}

inline void Document::read_string(const std::string& data,
                                  const std::string& name) {
  pegtl::parse_string<rules::file, Action, Errors>(data, name, *this);
  after_read(name);
}

inline void Document::read_memory(const char* data, const size_t size,
                                  const char* name) {
  pegtl::parse_memory<rules::file, Action, Errors>(data, size, name, *this);
  after_read(name);
}

inline void Document::read_cstream(std::FILE *f, const char* name,
                                   size_t maximum) {
  pegtl::parse_cstream<rules::file, Action, Errors>(f, name, maximum, *this);
  after_read(name);
}

inline void Document::read_istream(std::istream &is, const char* name,
                                   size_t maximum) {
  pegtl::parse_istream<rules::file, Action, Errors>(is, name, maximum, *this);
  after_read(name);
}

inline void Document::after_read(const std::string& name) {
  source = name;
  check_duplicates(*this);
}


inline bool check_file_syntax(const std::string& filename, std::string* msg) {
  pegtl::file_parser fp(filename);
  try {
#ifdef CIF_VALIDATE_SHOW_TRACE
    return fp.parse<rules::file, pegtl::nothing, pegtl::tracer>();
#else
    return fp.parse<rules::file, pegtl::nothing, Errors>();
#endif
  } catch (pegtl::parse_error& e) {
    if (msg)
      *msg = e.what();
    return false;
  }
}

// This will be moved to a different file
inline ValueType infer_valtype_of_string(const std::string& val) {
  assert(!val.empty());
  if (val == ".")
    return ValueType::Dot;
  else if (val == "?")
    return ValueType::QuestionMark;
  else if (is_numb(val))
    return ValueType::Numb;
  else
    return ValueType::Char;
}

inline void infer_valtypes_in_items(std::vector<Item>& items) {
  for (Item& item : items)
      if (item.type == ItemType::Value) {
        item.valtype = infer_valtype_of_string(item.tv.value);
      } else if (item.type == ItemType::Loop) {
        for (size_t i = 0; i != item.loop.tags.size(); ++i) {
          ValueType& vt = item.loop.tags[i].valtype;
          for (const std::string& v : LoopColumn{&item.loop, i}) {
            ValueType this_vt = infer_valtype_of_string(v);
            if (this_vt != vt) {
              // if we are here: vt != ValueType::Char
              if (vt == ValueType::NotSet || this_vt == ValueType::Numb) {
                vt = this_vt;
              } else if (this_vt == ValueType::Char) {
                vt = this_vt;
                break;
              }
            }
          }
        }
      } else if (item.type == ItemType::Frame) {
        infer_valtypes_in_items(item.frame.items);
      }
}

inline void infer_valtypes(Document &d) {
  for (Block& block : d.blocks)
    infer_valtypes_in_items(block.items);
}

inline bool is_text_field(const std::string& val) {
  size_t len = val.size();
  return len > 3 && val[0] == ';' && (val[len-2] == '\n' || val[len-2] == '\r');
}

} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
