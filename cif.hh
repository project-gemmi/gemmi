// Copyright 2017 Global Phasing Ltd.
#ifndef GEMMI_CIF_HH_
#define GEMMI_CIF_HH_
#include <cassert>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <new>
#include <unordered_set>

#include <pegtl.hh>
#ifdef CIF_VALIDATE_SHOW_TRACE
#include <pegtl/trace.hh>
#endif

namespace cif {

// the numeric type in CIF is called numb - it's a number with optional
// standard uncertainty in brackets: 1.23(8).
namespace numb_rules {
  using namespace pegtl;

  struct sign : opt<one<'+', '-'>> {};
  struct e : one<'e', 'E'> {};
  struct exponent : seq<sign, plus<digit>> {};
  struct uint_digit : digit {};
  struct fraction : plus<digit> {};
  struct base : if_then_else<one<'.'>, fraction,
                                       seq<plus<uint_digit>,
                                           opt<one<'.'>, opt<fraction>>>> {};
  // Error in brackets ,as per CIF spec. We ignore the value for now.
  struct err : seq<one<'('>, plus<digit>, one<')'>> {};
  struct numb : seq<sign, base, opt<e, exponent>, opt<err>, eof> {};
}

// Actions for getting the number. For now we ignore s.u., so the actions
// are equivalent to locale-independent std::stod().
template<typename Rule> struct ActionNumb : pegtl::nothing<Rule> {};
template<> struct ActionNumb<numb_rules::uint_digit> {
  template<typename Input> static void apply(const Input& in, double& d) {
      d = d * 10 + (*in.begin() - '0');
  }
};
template<> struct ActionNumb<numb_rules::fraction> {
  template<typename Input> static void apply(const Input& in, double& d) {
    double mult = 0.1;
    for (const auto* p = in.begin(); p != in.end(); ++p, mult *= 0.1)
      d += mult * (*p - '0');
  }
};
template<> struct ActionNumb<numb_rules::exponent> {
  template<typename Input> static void apply(const Input& in, double& d) {
    int n = 0;
    bool neg = false;
    const auto* p = in.begin();
    if (*p == '-')
      neg = true;
    else if (*p != '+')
      n = *p - '0';
    for (++p; p != in.end(); ++p)
      n = n * 10 + (*p - '0');
    // We don't expect too many exponents in CIF files, let's have
    // only this small LUT.
    static const double e[] = { 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };
    double uexp = n <= 8 ? e[n] : std::pow(10, n);
    if (neg)
      d /= uexp;
    else
      d *= uexp;
  }
};
template<> struct ActionNumb<numb_rules::numb> {
  template<typename Input> static void apply(const Input& in, double& d) {
    if (*in.begin() == '-')
      d = -d;
  }
};

inline bool is_numb(const std::string& s) {
  return pegtl::parse_string<numb_rules::numb, pegtl::nothing>(s, "");
}

inline double as_number(const std::string& s) {
  double d = 0;
  if (pegtl::parse_string<numb_rules::numb, ActionNumb>(s, "", d))
    return d;
  return NAN;
}



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
  struct str_data : pegtl_istring_t("data_") {};
  struct str_loop : pegtl_istring_t("loop_") {};
  struct str_global : pegtl_istring_t("global_") {};
  struct str_save : pegtl_istring_t("save_") {};
  struct str_stop : pegtl_istring_t("stop_") {};
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
  struct value: sor<singlequoted, doublequoted, textfield, unquoted> {};
  struct loop_tag : tag {};
  struct loop_value : value {};
  struct loop: if_must<str_loop, plus<seq<whitespace, loop_tag>>,
                                 star<seq<whitespace, loop_value>>,
                                 opt<seq<whitespace, str_stop>>> {};
  struct dataitem: if_must<tag, whitespace, value> {};
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
private:
  const std::string* cur_;
  const std::string* end_;
  int stride_;
};

struct LoopIter {
  const StrideIter begin_, end_;
  LoopIter() {}
  LoopIter(const Loop *loop, int idx)
    : begin_(loop->values, idx, loop->tags.size()),
      end_(loop->values, loop->values.size(), loop->tags.size()) {}
  StrideIter begin() const { return begin_; }
  StrideIter end() const { return end_; }
};


struct Item;

struct Block {
  std::string name;
  std::vector<Item> items;

  explicit Block(const std::string& name_) : name(name_) {}
  Block() {}

  const std::string* by_tag(const std::string& tag) const;
  const Loop* by_loop_tag(const std::string& tag, int* idx) const;
  LoopIter looped_values(const std::string& tag) const;
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


inline const std::string* Block::by_tag(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Value && i.tv.tag == tag)
      return &i.tv.value;
  return nullptr;
}

inline const Loop* Block::by_loop_tag(const std::string& tag, int* idx) const {
  for (const Item& i : items)
    if (i.type == ItemType::Loop) {
      int pos = i.loop.find_tag(tag);
      if (pos != -1) {
        if (idx)
          *idx = pos;
        return &i.loop;
      }
    }
  return nullptr;
}

inline LoopIter Block::looped_values(const std::string& tag) const {
  int idx = 0;
  const Loop* loop = by_loop_tag(tag, &idx);
  return loop ? LoopIter(loop, idx) : LoopIter();
}


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
          for (const std::string& v : LoopIter(&item.loop, i)) {
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

struct Document {
  Document() : items{nullptr} {}
  void parse_file(const std::string& filename);
  void infer_valtypes() {
    for (Block& block : blocks)
      infer_valtypes_in_items(block.items);
  }

  std::string source;
  std::vector<Block> blocks;
  std::vector<Comment> comments;
  std::vector<Item>* items; // items of the currently parsed block or frame
};

// **** parsing actions that fill the storage ****

template<typename Rule> struct Action : pegtl::nothing<Rule> {};

template<> struct Action<rules::datablockname> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.blocks.emplace_back(in.string());
    out.items = &out.blocks.back().items;
  }
};
template<> struct Action<rules::str_global> {
  template<typename Input> static void apply(const Input&, Document& out) {
    out.blocks.emplace_back();
    out.items = &out.blocks.back().items;
  }
};
template<> struct Action<rules::framename> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.items->emplace_back(in.string(), 0);
    out.items->back().line_number = in.line();
    out.items = &out.items->back().frame.items;
  }
};
template<> struct Action<rules::endframe> {
  template<typename Input> static void apply(const Input&, Document& out) {
    out.items = &out.blocks.back().items;
  }
};
template<> struct Action<rules::tag> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.items->emplace_back(in.string());
    out.items->back().line_number = in.line();
  }
};
template<> struct Action<rules::value> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items->back();
    assert(last_item.type == ItemType::Value);
    last_item.tv.value = in.string();
  }
};
template<> struct Action<rules::str_loop> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.items->emplace_back(0);
    out.items->back().line_number = in.line();
  }
};
template<> struct Action<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items->back();
    assert(last_item.type == ItemType::Loop);
    last_item.loop.tags.emplace_back(in.string());
    last_item.loop.tags.back().line_number = in.line();
  }
};
template<> struct Action<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items->back();
    assert(last_item.type == ItemType::Loop);
    last_item.loop.values.emplace_back(in.string());
  }
};
template<> struct Action<rules::loop> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    Item& last_item = out.items->back();
    assert(last_item.type == ItemType::Loop);
    const Loop& loop = last_item.loop;
    if (loop.values.size() % loop.tags.size() != 0)
      throw pegtl::parse_error("Wrong number of values in the loop", in);
  }
};
template<> struct Action<rules::comment> {
  template<typename Input> static void apply(const Input& in, Document& out) {
    out.comments.emplace_back(in.line(), in.string());
  }
};


void throw_validation_err(const Document& d, const Block& b, const Item& item,
                          const std::string& s) {
  throw std::runtime_error(d.source + ":" + std::to_string(item.line_number) +
                           " in data_" + b.name + ": " + s);
}

inline void check_no_name_dups(const Document& d) {
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

inline void Document::parse_file(const std::string& filename) {
  pegtl::parse_file<rules::file, Action, Errors>(filename, *this);
  source = filename;
  check_no_name_dups(*this);
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

} // namespace cif
#endif
// vim:sw=2:ts=2:et
