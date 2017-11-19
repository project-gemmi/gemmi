// Copyright 2017 Global Phasing Ltd.

#ifndef GEMMI_TO_CIF_HPP_
#define GEMMI_TO_CIF_HPP_

#include <fstream>
#include "cifdoc.hpp"

inline
std::ostream& operator<<(std::ostream& os, const gemmi::cif::Item& item) {
  using namespace gemmi::cif;
  switch (item.type) {
    case ItemType::Value:
      os << item.pair[0];
      os.put(is_text_field(item.pair[1]) ? '\n' : ' ');
      os << item.pair[1] << '\n';
      break;
    case ItemType::Loop: {
      if (item.loop.values.empty())
        break;
      os << "loop_";
      for (const LoopTag& tag : item.loop.tags)
        os << '\n' << tag.tag;
      size_t ncol = item.loop.tags.size();
      size_t col = 0;
      for (const std::string& val : item.loop.values) {
        os.put(col++ == 0 || is_text_field(val) ? '\n' : ' ');
        os << val;
        if (col == ncol)
          col = 0;
      }
      os << "\n\n"; // extra blank line for readability
      break;
    }
    case ItemType::Frame:
      os << "save_" << item.frame.name << '\n';
      for (const Item& inner_item : item.frame.items)
        os << inner_item;
      os << "save_\n";
      break;
    case ItemType::Comment:
      os << item.pair[1] << '\n';
      break;
    case ItemType::Erased:
      break;
  }
   return os;
}

inline
std::ostream& operator<<(std::ostream& os, const gemmi::cif::Document& doc) {
  bool first = true;
  for (const gemmi::cif::Block& block : doc.blocks) {
    if (!first)
      os.put('\n'); // extra blank line for readability
    os << "data_" << block.name << '\n';
    for (const gemmi::cif::Item& item : block.items)
      os << item;
    first = false;
  }
  return os;
}

namespace gemmi {
namespace cif {
inline void write_to_file(const Document& doc, const std::string& filename) {
  std::ofstream of(filename);
  if (!of)
    throw std::runtime_error("Failed to open " + filename);
  of << doc;
  of.close();
}
} // namespace cif
} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
