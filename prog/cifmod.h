
#include "gemmi/cifdoc.hpp"
#include <algorithm>           // for sort

enum CifModOptions { SkipCat=4, SortCif, AfterCifModOptions };

const option::Descriptor CifModUsage[] = {
  { 0, 0, 0, 0, 0, 0 }, // The first 3 entries are empty to make MapUsage[enum]
  { 0, 0, 0, 0, 0, 0 }, // work as expected.
  { 0, 0, 0, 0, 0, 0 }, // NoOp=0, Help=1, Version=2, Verbose=3.
  { 0, 0, 0, 0, 0, 0 },
  { SkipCat, 0, "", "skip-category", Arg::Required,
    "  --skip-category=CAT  \tDo not output tags starting with _CAT" },
  { SortCif, 0, "", "sort", Arg::None,
    "  --sort  \tSort tags in alphabetical order." },
};

inline const std::string& sort_name(const gemmi::cif::Item& item) {
  static const std::string none;
  switch (item.type) {
    case gemmi::cif::ItemType::Pair:
      return item.pair[0];
    case gemmi::cif::ItemType::Loop:
      return item.loop.tags.empty() ? none : item.loop.tags[0];
    case gemmi::cif::ItemType::Frame:
      return item.frame.name;
    case gemmi::cif::ItemType::Comment:
    case gemmi::cif::ItemType::Erased:
      return none;
  }
  gemmi::unreachable();
}

inline void sort_items(gemmi::cif::Block& block) {
  std::sort(block.items.begin(), block.items.end(),
            [&](const gemmi::cif::Item& a, const gemmi::cif::Item& b) {
                return sort_name(a) < sort_name(b);
            });
  for (gemmi::cif::Item& item : block.items)
    if (item.type == gemmi::cif::ItemType::Loop) {
      std::vector<std::string>& tags = item.loop.tags;
      std::vector<std::string>& values = item.loop.values;
      int n = (int) tags.size();
      std::vector<int> new_index(n);
      for (int i = 0; i != n; ++i)
        new_index[i] = i;
      std::sort(new_index.begin(), new_index.end(), [&](int a, int b) {
          return gemmi::to_lower(tags[a]) < gemmi::to_lower(tags[b]);
      });
      std::vector<std::string> tmp = tags;
      for (int i = 0; i != n; ++i)
        tags[i] = tmp[new_index[i]];
      for (size_t offset = 0; offset != values.size(); offset += n) {
        for (int i = 0; i != n; ++i)
          tmp[i] = values[offset + i];
        for (int i = 0; i != n; ++i)
          values[offset + i] = tmp[new_index[i]];
      }
    }
}

inline void apply_cif_doc_modifications(gemmi::cif::Document& doc,
                                   const std::vector<option::Option>& options) {
  for (const option::Option* opt = options[SkipCat]; opt; opt = opt->next()) {
    std::string category = gemmi::to_lower(opt->arg);
    if (category[0] != '_')
      category.insert(0, 1, '_');
    for (gemmi::cif::Block& block : doc.blocks)
      for (gemmi::cif::Item& item : block.items)
        if (item.has_prefix(category))
          item.erase();
  }
  if (options[SortCif])
    for (gemmi::cif::Block& block : doc.blocks)
      sort_items(block);
}

