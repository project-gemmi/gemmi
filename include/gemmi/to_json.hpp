// Copyright 2017 Global Phasing Ltd.

#ifndef GEMMI_TO_JSON_HPP_
#define GEMMI_TO_JSON_HPP_
#include "cif.hpp"
#include <cassert>
#include <ostream>

namespace gemmi {
namespace cif {

class JsonWriter {
public:
  bool use_bare_tags = false;  // "tag" instead of "_tag"
  int quote_numbers = 1;  // 0=never (no s.u.), 1=mix, 2=always
  std::string unknown = "null";  // how to convert '?' from CIF
  explicit JsonWriter(std::ostream& os) : os_(os), linesep_("\n ") {}
  void write_json(const Document& d);

private:
  std::ostream& os_;
  std::string linesep_;

  // based on tao/json/internal/escape.hpp
  static void escape(std::ostream& os, const std::string& s, size_t pos,
                     bool to_lower) {
    static const char* h = "0123456789abcdef";
    const char* p = s.data() + pos;
    const char* l = p;
    const char* const e = s.data() + s.size();
    while (p != e) {
      const unsigned char c = *p;
      if (c == '\\') {
        os.write(l, p - l);
        l = ++p;
        os << "\\\\";
      } else if (c == '"') {
        os.write(l, p - l);
        l = ++p;
        os << "\\\"";
      } else if (c < 32) {
        os.write(l, p - l);
        l = ++p;
        switch ( c ) {
          case '\b': os << "\\b"; break;
          case '\f': os << "\\f"; break;
          case '\n': os << "\\n"; break;
          case '\r': os << "\\r"; break;
          case '\t': os << "\\t"; break;
          default: os << "\\u00" << h[(c & 0xf0) >> 4] << h[c & 0x0f];
        }
      } else if (to_lower && c >= 'A' && c <= 'Z') {
        os.write(l, p - l);
        l = ++p;
        os.put(c + 32);
      } else if (c == 127) {
        os.write(l, p - l);
        l = ++p;
        os << "\\u007f";
      } else {
        ++p;
      }
    }
    os.write(l, p - l);
  }

  void write_string(const std::string& s, size_t pos=0, bool to_lower=false) {
    os_.put('"');
    escape(os_, s, pos, to_lower);
    os_.put('"');
  }

  void write_tag(const std::string& tag) {
    write_string(tag, use_bare_tags ? 1 : 0, true);
  }

  void write_as_number(const std::string& value) {
    // if we are here, value is not empty
    if (value[0] == '.') // in JSON numbers cannot start with dot
      os_.put('0');
    // in JSON the number cannot start with +
    size_t pos = 0;
    if (value[pos] == '+') {
      pos = 1;
    } else if (value[pos] == '-') { // make handling -001 easier
      os_.put('-');
      pos = 1;
    }
    // in JSON left-padding with 0s is not allowed
    while (value[pos] == '0' && isdigit(value[pos+1]))
      ++pos;
    // in JSON dot must be followed by digit
    size_t dotpos = value.find('.');
    if (dotpos != std::string::npos && !isdigit(value[dotpos+1])) {
      os_ << value.substr(pos, dotpos+1-pos) << '0';
      pos = dotpos + 1;
    }
    if (value.back() != ')')
      os_ << value.c_str() + pos;
    else
      os_ << value.substr(pos, value.find('(', pos) - pos);
  }

  void write_value(const std::string& value) {
    if (value == "?")
      os_ << unknown;
    else if (value == ".")
      os_ << "null";
    else if (quote_numbers < 2 && is_numb(value) &&
             (quote_numbers == 0 || value.back() != ')'))
      write_as_number(value);
    else
      write_string(as_string(value));
  }

  void write_item(const Item& item) {
    switch (item.type) {
      case ItemType::Value:
        write_tag(item.tv.tag);
        os_ << ": ";
        write_value(item.tv.value);
        break;
      case ItemType::Loop: {
        size_t ncol = item.loop.tags.size();
        const auto& vals = item.loop.values;
        for (size_t i = 0; i < ncol; i++) {
          if (i != 0)
            os_ << "," << linesep_;
          write_tag(item.loop.tags[i].tag);
          os_ << ": [";
          for (size_t j = i; j < vals.size(); j += ncol) {
            if (j != i)
              os_.put(',');
            write_value(vals[j]);
          }
          os_.put(']');
        }
        break;
      }
      case ItemType::Frame:
        write_map(item.frame.name, item.frame.items);
        break;
      case ItemType::Comment:
        break;
      case ItemType::Erased:
        break;
    }
  }

  // works for both block and frame
  void write_map(const std::string& name, const std::vector<Item>& items) {
    write_string(name, 0, true);
    os_ << ": ";
    size_t n = linesep_.size();
    linesep_.resize(n + 1, ' ');
    char first = '{';
    for (const Item& item : items) {
      os_ << first << linesep_;
      write_item(item);
      first = ',';
    }
    linesep_.resize(n);
    os_ << linesep_ << '}';
  }
};

inline void JsonWriter::write_json(const Document& d) {
  os_.put('{');
  for (const Block& block : d.blocks) {
    if (&block != &d.blocks[0])
      os_.put(',');
    os_ << linesep_;
    write_map(block.name, block.items);
  }
  os_ << "\n}\n";
}

} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
