// Copyright 2017 Global Phasing Ltd.

#ifndef GEMMI_TO_CIF_HPP_
#define GEMMI_TO_CIF_HPP_

#include <fstream>
#include "cifdoc.hpp"

namespace gemmi {
namespace cif {

enum class Style {
  Simple,
  NoBlankLines,
  PreferPairs,  // write single-row loops as pairs
  Pdbx,         // PreferPairs + put '#' (empty comments) between categories
};

inline void write_out_pair(std::ostream& os,
                           const std::string& name, const std::string& value) {
  bool br = is_text_field(value) || name.size() + value.size() > 120;
  os << name << (br ? '\n' : ' ') << value << '\n';
}

inline void write_out_loop(std::ostream& os, const Loop& loop, Style style) {
  if (loop.values.empty())
    return;
  if ((style == Style::PreferPairs || style == Style::Pdbx) &&
      loop.length() == 1) {
    for (size_t i = 0; i != loop.tags.size(); ++i)
      write_out_pair(os, loop.tags[i], loop.values[i]);
    return;
  }
  os << "loop_";
  for (const std::string& tag : loop.tags)
    os << '\n' << tag;
  size_t ncol = loop.tags.size();
  size_t col = 0;
  for (const std::string& val : loop.values) {
    os.put(col++ == 0 || is_text_field(val) ? '\n' : ' ');
    os << val;
    if (col == ncol)
      col = 0;
  }
  os.put('\n');
}

inline void write_out_item(std::ostream& os, const Item& item, Style style) {
  switch (item.type) {
    case ItemType::Pair:
      write_out_pair(os, item.pair[0], item.pair[1]);
      break;
    case ItemType::Loop:
      write_out_loop(os, item.loop, style);
      break;
    case ItemType::Frame:
      os << "save_" << item.frame.name << '\n';
      for (const Item& inner_item : item.frame.items)
        write_out_item(os, inner_item, style);
      os << "save_\n";
      break;
    case ItemType::Comment:
      os << item.pair[1] << '\n';
      break;
    case ItemType::Erased:
      break;
  }
}

inline bool should_be_separted_(const Item& a, const Item& b) {
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

inline void write_cif_to_stream(std::ostream& os, const Document& doc,
                                Style s) {
  bool first = true;
  for (const Block& block : doc.blocks) {
    if (!first)
      os.put('\n'); // extra blank line for readability
    os << "data_" << block.name << '\n';
    if (s == Style::Pdbx)
      os << "#\n";
    const Item* prev = nullptr;
    for (const Item& item : block.items)
      if (item.type != ItemType::Erased) {
        if (prev && s != Style::NoBlankLines &&
            should_be_separted_(*prev, item))
          os << (s == Style::Pdbx ? "#\n" : "\n");
        write_out_item(os, item, s);
        prev = &item;
      }
    first = false;
    if (s == Style::Pdbx)
      os << "#\n";
  }
}

inline void write_cif_to_file(const Document& doc, const std::string& filename,
                          Style s=Style::Simple) {
  std::ofstream of(filename);
  if (!of)
    throw std::runtime_error("Failed to open " + filename);
  write_cif_to_stream(of, doc, s);
  of.close();
}

} // namespace cif
} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
