// Copyright 2017 Global Phasing Ltd.

// Writing cif::Document or its parts to std::ostream.

#ifndef GEMMI_TO_CIF_HPP_
#define GEMMI_TO_CIF_HPP_

#include <ostream>
#include "cifdoc.hpp"

namespace gemmi {
namespace cif {

enum class Style {
  Simple,
  NoBlankLines,
  PreferPairs,  // write single-row loops as pairs
  Pdbx,         // PreferPairs + put '#' (empty comments) between categories
  Indent35,     // start values in pairs from 35th column
  Aligned,      // columns in tables are left-aligned
};

struct WriteOptions {
  /// write single-row loops as pairs
  bool prefer_pairs = false;
  /// no blank lines between categories, only between blocks
  bool compact = false;
  /// put '#' (empty comments) before/after categories
  bool misuse_hash = false;
  /// width reserved for tags in pairs (e.g. 34 = value starts at 35th column)
  std::uint16_t align_pairs = 0;
  /// if non-zero, determines max width of each column in a loop and aligns
  /// all values to this width; the width is capped with the given value
  std::uint16_t align_loops = 0;

  WriteOptions() {}
  // implicit conversion from deprecated Style (for backward compatibility)
  WriteOptions(Style style) {
    switch (style) {
      case Style::Simple:
        break;
      case Style::NoBlankLines:
        compact = true;
        break;
      case Style::PreferPairs:
        prefer_pairs = true;
        break;
      case Style::Pdbx:
        prefer_pairs = true;
        misuse_hash = true;
        break;
      case Style::Indent35:
        align_pairs = 33;
        break;
      case Style::Aligned:
        align_pairs = 33;
        align_loops = 30;
        break;
    }
  }
  std::string str() const {
    std::string s;
    if (prefer_pairs)
      s += "prefer_pairs,";
    if (compact)
      s += "compact,";
    if (misuse_hash)
      s += "misuse_hash,";
    if (align_pairs != 0)
      s += "align_pairs=" + std::to_string(align_pairs) + ",";
    if (align_loops != 0)
      s += "align_loops=" + std::to_string(align_loops) + ",";
    if (!s.empty())
      s.pop_back();
    return s;
  }
};

/// std::ostream with buffering. C++ streams are so slow that even primitive
/// buffering makes it significantly more efficient.
class BufOstream {
public:
  explicit BufOstream(std::ostream& os_) : os(os_), ptr(buf) {}
  ~BufOstream() { flush(); }
  void flush() {
    os.write(buf, ptr - buf);
    ptr = buf;
  }
  void write(const char* s, size_t len) {
    constexpr int margin = sizeof(buf) - 512;
    if (ptr - buf + len > margin) {
      flush();
      if (len > margin) {
        os.write(s, len);
        return;
      }
    }
    std::memcpy(ptr, s, len);
    ptr += len;
  }
  void operator<<(const std::string& s) {
    write(s.c_str(), s.size());
  }
  // below we don't check the buffer boundary, these functions add <512 bytes
  void put(char c) {
    *ptr++ = c;
  }
  void pad(size_t n) {
    std::memset(ptr, ' ', n);
    ptr += n;
  }
private:
  std::ostream& os;
  // increasing buffer to 8kb or 64kb doesn't make significant difference
  char buf[4096];
  char* ptr;
};

// CIF files are read in binary mode. It makes difference only for text fields.
// If the text field with \r\n would be written as is in text mode on Windows
// \r would get duplicated. As a workaround, here we convert \r\n to \n.
// Hopefully \r that gets removed here is never meaningful.
inline void write_text_field(BufOstream& os, const std::string& value) {
  for (size_t pos = 0, end = 0; end != std::string::npos; pos = end + 1) {
    end = value.find("\r\n", pos);
    size_t len = (end == std::string::npos ? value.size() : end) - pos;
    os.write(value.c_str() + pos, len);
  }
}

inline void write_out_pair(BufOstream& os, const std::string& name,
                           const std::string& value, WriteOptions options) {
  os << name;
  if (is_text_field(value)) {
    os.put('\n');
    write_text_field(os, value);
  } else {
    if (name.size() + value.size() > 120) {
      os.put('\n');
    } else {
      os.put(' ');
      if (name.size() < options.align_pairs)
        os.pad(options.align_pairs - name.size());
    }
    os << value;
  }
  os.put('\n');
}

inline void write_out_loop(BufOstream& os, const Loop& loop, WriteOptions options) {
  if (loop.values.empty())
    return;
  if (options.prefer_pairs && loop.length() == 1) {
    for (size_t i = 0; i != loop.tags.size(); ++i)
      write_out_pair(os, loop.tags[i], loop.values[i], options);
    return;
  }
  // tags
  os.write("loop_", 5);
  for (const std::string& tag : loop.tags) {
    os.put('\n');
    os << tag;
  }
  // values
  size_t ncol = loop.tags.size();

  std::vector<size_t> col_width(ncol, 0);
  if (options.align_loops > 0) {
    size_t col = 0;
    for (const std::string& val : loop.values) {
      if (!is_text_field(val))
        col_width[col] = std::max(col_width[col], val.size());
      if (++col == ncol)
        col = 0;
    }
    for (size_t& w : col_width)
      w = std::min(w, (size_t)options.align_loops);
  }

  size_t col = 0;
  bool need_new_line = true;
  for (const std::string& val : loop.values) {
    bool text_field = is_text_field(val);
    os.put(need_new_line || text_field ? '\n' : ' ');
    need_new_line = text_field;
    if (text_field)
      write_text_field(os, val);
    else
      os << val;
    if (col != ncol - 1) {
      if (val.size() < col_width[col])
        os.pad(col_width[col] - val.size());
      ++col;
    } else {
      col = 0;
      need_new_line = true;
    }
  }
  os.put('\n');
}

inline void write_out_item(BufOstream& os, const Item& item, WriteOptions options) {
  switch (item.type) {
    case ItemType::Pair:
      write_out_pair(os, item.pair[0], item.pair[1], options);
      break;
    case ItemType::Loop:
      write_out_loop(os, item.loop, options);
      break;
    case ItemType::Frame:
      os.write("save_", 5);
      os << item.frame.name;
      os.put('\n');
      for (const Item& inner_item : item.frame.items)
        write_out_item(os, inner_item, options);
      os.write("save_\n", 6);
      break;
    case ItemType::Comment:
      os << item.pair[1];
      os.put('\n');
      break;
    case ItemType::Erased:
      break;
  }
}

inline bool should_be_separated_(const Item& a, const Item& b) {
  if (a.type == ItemType::Comment || b.type == ItemType::Comment)
    return false;
  if (a.type != ItemType::Pair || b.type != ItemType::Pair)
    return true;
  // check if we have mmcif-like tags from different categories
  auto adot = a.pair[0].find('.');
  if (adot == std::string::npos)
    return false;
  auto bdot = b.pair[0].find('.');
  return adot != bdot || a.pair[0].compare(0, adot, b.pair[0], 0, adot) != 0;
}

inline void write_cif_block_to_stream(std::ostream& os_, const Block& block,
                                      WriteOptions options=WriteOptions()) {
  BufOstream os(os_);
  os.write("data_", 5);
  os << block.name;
  os.put('\n');
  if (options.misuse_hash)
    os.write("#\n", 2);
  const Item* prev = nullptr;
  for (const Item& item : block.items) {
    if (item.type == ItemType::Erased)
      continue;
    if (prev && !options.compact && should_be_separated_(*prev, item)) {
      if (options.misuse_hash)
        os.put('#');
      os.put('\n');
    }
    write_out_item(os, item, options);
    prev = &item;
  }
  if (options.misuse_hash)
    os.write("#\n", 2);
}

inline void write_cif_to_stream(std::ostream& os, const Document& doc,
                                WriteOptions options=WriteOptions()) {
  bool first = true;
  for (const Block& block : doc.blocks) {
    if (!first)
      os.put('\n'); // extra blank line for readability
    write_cif_block_to_stream(os, block, options);
    first = false;
  }
}

} // namespace cif
} // namespace gemmi

#endif
