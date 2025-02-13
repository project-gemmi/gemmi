// Copyright 2017 Global Phasing Ltd.

#include <sstream>
#include "gemmi/cif.hpp"  // for parse_input
#include "gemmi/cifdoc.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/fstream.hpp"
#include "gemmi/ddl.hpp"
#include "gemmi/read_cif.hpp"  // for read_cif_gz

#include "common.h"
#include "serial.h"  // for getstate, setstate
#include "make_iterator.h"
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

using namespace gemmi::cif;

namespace {

std::string pyobject_to_string(nb::handle handle, bool raw) {
  PyObject* ptr = handle.ptr();
  if (ptr == Py_None)
    return "?";
  if (ptr == Py_False)
    return ".";
  if (ptr == Py_True)
    throw nb::value_error("unexpected value True");
  if (raw || PyFloat_Check(ptr) || PyLong_Check(ptr))
    return nb::str(handle).c_str();
  return gemmi::cif::quote(nb::str(handle).c_str());
}

nb::object string_to_pyobject(const std::string& s, bool raw) {
  if (raw)
    return nb::cast(s);
  if (is_null(s))
    return nb::borrow(s[0] == '?' ? Py_None : Py_False);
  return nb::cast(as_string(s));
}

std::vector<std::string> quote_list(const nb::list& items) {
  size_t size = items.size();
  std::vector<std::string> ret;
  ret.reserve(size);
  for (auto handle : items)
    ret.push_back(pyobject_to_string(handle, false));
  return ret;
}

template<typename T>
T& add_to_vector(std::vector<T>& vec, const T& new_item, int pos) {
  if (pos < 0)
    pos = (int) vec.size();
  else if (pos > (int) vec.size())
    throw nb::index_error();
  vec.insert(vec.begin() + pos, new_item);
  return vec[pos];
}

// cif parsing without checking for missing values and duplicates
void cif_parse_string(Document& doc, const std::string& data) {
  tao::pegtl::memory_input<> in(data, "string");
  parse_input(doc, in);
}
void cif_parse_file(Document& doc, const std::string& filename) {
  GEMMI_CIF_FILE_INPUT(in, filename);
  parse_input(doc, in);
}

}  // anonymous namespace

// for delitem_slice
template<>
void delitem_at_index(Table& t, size_t idx) { t.remove_row((int)idx); }
template<>
void delitem_range(Table& t, size_t start, size_t end) { t.remove_rows((int)start, (int)end); }

void add_cif(nb::module_& cif) {
  nb::class_<Block> cif_block(cif, "Block");
  nb::class_<Item> cif_item(cif, "Item");
  nb::class_<Loop> cif_loop(cif, "Loop");
  nb::class_<Column> cif_column(cif, "Column");
  nb::class_<Table> cif_table(cif, "Table");
  nb::class_<Table::Row> cif_table_row(cif_table, "Row");

  nb::enum_<Style>(cif, "Style")
    .value("Simple", Style::Simple)
    .value("NoBlankLines", Style::NoBlankLines)
    .value("PreferPairs", Style::PreferPairs)
    .value("Pdbx", Style::Pdbx)
    .value("Indent35", Style::Indent35)
    .value("Aligned", Style::Aligned)
    ;
  nb::class_<WriteOptions>(cif, "WriteOptions")
    .def(nb::init<>())
    .def(nb::init_implicit<Style>())
    .def_rw("prefer_pairs", &WriteOptions::prefer_pairs)
    .def_rw("compact", &WriteOptions::compact)
    .def_rw("misuse_hash", &WriteOptions::misuse_hash)
    .def_rw("align_pairs", &WriteOptions::align_pairs)
    .def_rw("align_loops", &WriteOptions::align_loops)
    .def("__repr__", [](const WriteOptions &self) -> std::string {
        std::string str = self.str();
        if (str.empty())
          return "gemmi.cif.WriteOptions()";
        return gemmi::cat("<gemmi.cif.WriteOptions ", str, '>');
    });
  nb::class_<Document>(cif, "Document")
    .def(nb::init<>())
    .def_rw("source", &Document::source)
    .def("__len__", [](const Document& d) { return d.blocks.size(); })
    .def("__iter__", [](Document& d) { return usual_iterator(d, d.blocks); },
         nb::keep_alive<0, 1>())
    .def("__getitem__", [](Document& d, const std::string& name) -> Block& {
        Block* b = d.find_block(name);
        if (!b)
          throw nb::key_error(gemmi::cat("block '", name, "' does not exist").c_str());
        return *b;
    }, nb::arg("name"), nb::rv_policy::reference_internal)
    .def("__getitem__", [](Document& d, int index) -> Block& {
        return d.blocks[normalize_index(index, d.blocks)];
    }, nb::arg("index"), nb::rv_policy::reference_internal)
    .def("__getitem__", [](Document &d, const nb::slice& slice) -> nb::list {
        return getitem_slice(d.blocks, slice);
    }, nb::rv_policy::reference_internal)
    .def("__contains__", [](Document& d, const std::string& name) -> bool {
        Block* b = d.find_block(name);
        return b != nullptr;
    }, nb::arg("name"))
    .def("__delitem__", [](Document &d, int index) {
        if (index < 0)
          index += (int) d.blocks.size();
        if (index < 0 || static_cast<size_t>(index) >= d.blocks.size())
          throw nb::index_error();
        d.blocks.erase(d.blocks.begin() + index);
    }, nb::arg("index"))
    .def("add_copied_block", [](Document& d, const Block& block, int pos) -> Block& {
        return add_to_vector(d.blocks, block, pos);
    }, nb::arg("block"), nb::arg("pos")=-1,
       nb::rv_policy::reference_internal)
    .def("add_new_block", &Document::add_new_block,
         nb::arg("name"), nb::arg("pos")=-1,
         nb::rv_policy::reference_internal)
    .def("clear", &Document::clear)
    .def("parse_string", &cif_parse_string)
    .def("parse_file", &cif_parse_file)
    .def("check_for_missing_values", &check_for_missing_values)
    .def("check_for_duplicates", &check_for_duplicates)
    .def("sole_block", (Block& (Document::*)()) &Document::sole_block,
         nb::rv_policy::reference_internal,
         "Returns the only block if there is exactly one")
    .def("find_block",
         (Block* (Document::*)(const std::string&)) &Document::find_block,
         nb::arg("name"),
         nb::rv_policy::reference_internal)
    .def("write_file",
         [](const Document& doc, const std::string& filename, WriteOptions opt) {
        gemmi::Ofstream os(filename);
        write_cif_to_stream(os.ref(), doc, opt);
    }, nb::arg("filename"), nb::arg("options")=WriteOptions(),
    "Write data to a CIF file.")
    .def("as_string", [](const Document& d, WriteOptions opt) {
        std::ostringstream os;
        write_cif_to_stream(os, d, opt);
        return os.str();
    }, nb::arg("options")=WriteOptions(), "Returns a string in CIF format.")
    .def("as_json", [](const Document& d, bool mmjson, bool lowercase_names) {
        std::ostringstream os;
        JsonWriteOptions options;
        if (mmjson) {
          options = JsonWriteOptions::mmjson();
        } else {
          // in C++17 std::optional<bool> would be used
          options.lowercase_names = lowercase_names;
        }
        write_json_to_stream(os, d, options);
        return os.str();
    }, nb::arg("mmjson")=false, nb::arg("lowercase_names")=true,
    "Returns JSON representation in a string.")
    .def("__getstate__", &getstate<Document>)
    .def("__setstate__", &setstate<Document>)
    .def("__repr__", [](const Document &d) {
        std::string s = "<gemmi.cif.Document with ";
        s += std::to_string(d.blocks.size());
        s += " blocks (";
        for (size_t i = 0; i != std::min(size_t{3}, d.blocks.size()); ++i) {
          if (i != 0)
            s += ", ";
          s += d.blocks[i].name;
        }
        s += d.blocks.size() > 3 ? "...)>" : ")>";
        return s;
    });

  cif_block
    .def(nb::init<const std::string &>())
    .def_rw("name", &Block::name)
    .def("__iter__", [](Block& self) { return usual_iterator(self, self.items); },
         nb::keep_alive<0, 1>())
    .def("__getitem__", [](Block& self, int index) -> Item& {
        return self.items[normalize_index(index, self.items)];
    }, nb::arg("index"), nb::rv_policy::reference_internal)
    .def("find_pair", [](const Block& self, const std::string& tag) -> nb::object {
        if (const Pair* p = self.find_pair(tag))
          return nb::make_tuple((*p)[0], (*p)[1]);
        return nb::none();
     }, nb::arg("tag"))
    .def("find_pair_item", &Block::find_pair_item, nb::arg("tag"),
         nb::rv_policy::reference_internal)
    .def("find_value", &Block::find_value, nb::arg("tag"),
         nb::rv_policy::reference)
    .def("find_loop", &Block::find_loop, nb::arg("tag"),
         nb::keep_alive<0, 1>())
    .def("find_loop_item", &Block::find_loop_item, nb::arg("tag"),
         nb::rv_policy::reference_internal)
    .def("find_values", &Block::find_values, nb::arg("tag"),
         nb::keep_alive<0, 1>())
    .def("find", (Table (Block::*)(const std::string&,
            const std::vector<std::string>&)) &Block::find,
         nb::arg("prefix"), nb::arg("tags"), nb::keep_alive<0, 1>())
    .def("find", (Table (Block::*)(const std::vector<std::string>&))
                 &Block::find,
         nb::arg("tags"), nb::keep_alive<0, 1>())
    .def("find_or_add", &Block::find_or_add,
         nb::arg("prefix"), nb::arg("tags"), nb::keep_alive<0, 1>())
    .def("add_item", [](Block& block, const Item& item, int pos) {
        return add_to_vector(block.items, item, pos);
    }, nb::arg("item"), nb::arg("pos")=-1,
       nb::rv_policy::reference_internal)
    .def("find_frame", &Block::find_frame, nb::arg("name"),
         nb::rv_policy::reference_internal)
    .def("item_as_table", &Block::item_as_table)
    .def("get_index", &Block::get_index, nb::arg("tag"))
    .def("set_pair", &Block::set_pair, nb::arg("tag"), nb::arg("value"))
    .def("set_pairs",
         [](Block &self, const std::string& prefix, const nb::dict& data, bool raw) {
           ItemSpan span(self.items, prefix);
           for (const auto& item : data) {
             std::string key = nb::str(item.first).c_str();
             std::string value = pyobject_to_string(item.second, raw);
             span.set_pair(prefix + key, value);
           }
         }, nb::arg("prefix"), nb::arg("data"), nb::arg("raw")=false)
    .def("init_loop", &Block::init_loop, nb::arg("prefix"), nb::arg("tags"),
         nb::rv_policy::reference_internal)
    .def("move_item", &Block::move_item, nb::arg("old_pos"), nb::arg("new_pos"))
    .def("find_mmcif_category", &Block::find_mmcif_category,
         nb::arg("category"), nb::keep_alive<0, 1>(),
         "Returns Table with all items in the category.")
    .def("get_mmcif_category_names", &Block::get_mmcif_category_names,
         "For mmCIF files only. Returns list of all category prefixes (_x.)")
    .def("init_mmcif_loop", &Block::init_mmcif_loop,
         nb::arg("cat"), nb::arg("tags"),
         nb::rv_policy::reference_internal)
    .def("set_mmcif_category",
         [](Block &self, std::string name, const nb::dict& data, bool raw) {
           size_t w = data.size();
           std::vector<std::string> tags;
           tags.reserve(w);
           std::vector<nb::list> values;
           values.reserve(w);
           for (const auto& item : data) {
             tags.emplace_back(nb::str(item.first).c_str());
             // We could support all iterables as values, but it'd be confusing
             // if strings/bytes were converted to sequences of letters.
             nb::list& ref = values.emplace_back(nb::cast<nb::list>(item.second));
             if (ref.size() != values[0].size())
               throw nb::value_error("all columns must have equal length");
           }
           Loop& loop = self.init_mmcif_loop(std::move(name), std::move(tags));
           loop.values.resize(w != 0 ? w * values[0].size() : 0);
           for (size_t col = 0; col != w; ++col) {
             size_t idx = col;
             for (auto handle : values[col]) {
               loop.values[idx] = pyobject_to_string(handle, raw);
               idx += w;
             }
           }
         }, nb::arg("name"), nb::arg("data"), nb::arg("raw")=false)
    .def("get_mmcif_category",
         [](Block &self, std::string name, bool raw) {
           ensure_mmcif_category(name);
           nb::dict data;
           Table table = self.find_mmcif_category(name);
           int len = (int) table.length();
           for (const std::string& tag : table.tags()) {
             assert(tag.size() >= name.size());
             Column col = table.column((int)data.size());
             // We could avoid reallocations with:
             //  nb::list new_list = nb::steal<nb::list>(PyList_New(len));
             // but then would need to use directly NB_LIST_SET_ITEM
             // because new_list[i]=... doesn't check for null.
             nb::list new_list;
             for (int i = 0; i != len; ++i)
               new_list.append(string_to_pyobject(col[i], raw));
             data[tag.c_str() + name.size()] = new_list;
           }
           return data;
         }, nb::arg("name"), nb::arg("raw")=false)
    .def("write_file",
         [](const Block& self, const std::string& filename, WriteOptions opt) {
        gemmi::Ofstream os(filename);
        write_cif_block_to_stream(os.ref(), self, opt);
    }, nb::arg("filename"), nb::arg("options")=WriteOptions(),
    "Write data to a CIF file.")
    .def("as_string", [](const Block& self, WriteOptions opt) {
        std::ostringstream os;
        write_cif_block_to_stream(os, self, opt);
        return os.str();
    }, nb::arg("options")=WriteOptions(), "Returns a string in CIF format.")
    .def("__getstate__", &getstate<Block>)
    .def("__setstate__", &setstate<Block>)
    .def("__repr__", [](const Block &self) {
        return gemmi::cat("<gemmi.cif.Block ", self.name, '>');
    });


  cif_item
    .def("erase", &Item::erase)
    .def_ro("line_number", &Item::line_number)
    .def_prop_ro("pair", [](Item& self) -> nb::object {
        if (self.type == ItemType::Pair)
          return nb::make_tuple(self.pair[0], self.pair[1]);
        return nb::none();
    })
    .def_prop_ro("loop", [](Item& self) {
        return self.type == ItemType::Loop ? &self.loop : nullptr;
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("frame", [](Item& self) {
        return self.type == ItemType::Frame ? &self.frame : nullptr;
    }, nb::rv_policy::reference_internal)
    ;

  cif_loop
    .def(nb::init<>())
    .def("width", &Loop::width, "Returns number of columns")
    .def("length", &Loop::length, "Returns number of rows")
    .def_ro("tags", &Loop::tags)
    .def_ro("values", &Loop::values)
    .def("__getitem__", [](Loop &self, std::pair<int,int> idx) {
        return self.val(c_index(idx.first, self.length()),
                        c_index(idx.second, self.width()));
    })
    .def("__setitem__", [](Loop &self, std::pair<int,int> idx, const std::string& value) {
        self.val(c_index(idx.first, self.length()),
                 c_index(idx.second, self.width())) = value;
    })
    .def("add_row", &Loop::add_row<std::vector<std::string>>,
         nb::arg("new_values"), nb::arg("pos")=-1)
    .def("add_columns", &Loop::add_columns,
         nb::arg("column_names"), nb::arg("value"), nb::arg("pos")=-1)
    .def("remove_column", &Loop::remove_column)
    .def("set_all_values", &Loop::set_all_values, nb::arg("columns"))
    .def("__repr__", [](const Loop &self) {
        return gemmi::cat("<gemmi.cif.Loop ", self.length(), " x ", self.width(), '>');
    });


  cif_column
    .def(nb::init<>())
    .def("get_loop", &Column::get_loop,
         nb::rv_policy::reference_internal)
    .def_prop_rw("tag",
                  [](const Column& self) { return self.get_tag(); },
                  [](Column& self, const std::string& s) {
                      if (std::string* tag = self.get_tag())
                          *tag = s;
                  })
    .def("__iter__", [](Column& self) { return usual_iterator(self, self); },
         nb::keep_alive<0, 1>())
    .def("__bool__", [](const Column &self) { return self.item() != nullptr; })
    .def("__len__", [](const Column &self) { return self.length(); })
    .def("__getitem__", (std::string& (Column::*)(int)) &Column::at)
    .def("__setitem__", [](Column &self, int idx, const std::string& value) {
        self.at(idx) = value;
    })
    .def("erase", &Column::erase)
    .def("__repr__", [](const Column &self) {
        std::string s = "<gemmi.cif.Column ";
        if (const std::string* tag = self.get_tag())
          gemmi::cat_to(s, *tag, " length ", self.length(), '>');
        else
          s += "nil>";
        return s;
    })
    .def("str", &Column::str, nb::arg("index"));

  cif_table
    .def("width", &Table::width)
    .def_ro("prefix_length", &Table::prefix_length)
    .def_prop_ro("loop", &Table::get_loop,
                           nb::rv_policy::reference_internal)
    .def("get_prefix", &Table::get_prefix)
    .def("has_column", &Table::has_column)
    .def("column", &Table::column, nb::arg("n"), nb::keep_alive<0, 1>())
    .def("find_row", &Table::find_row, nb::keep_alive<0, 1>())
    .def("find_column", &Table::find_column, nb::arg("tag"),
         nb::keep_alive<0, 1>())
    .def("erase", &Table::erase)
    .def("ensure_loop", &Table::ensure_loop)
    .def("append_row", &Table::append_row<std::vector<std::string>>,
         nb::arg("new_values"))
    .def("remove_row", &Table::remove_row, nb::arg("row_index"))
    .def("move_row", &Table::move_row, nb::arg("old_pos"), nb::arg("new_pos"))
    .def_prop_ro("tags", &Table::tags, nb::keep_alive<0, 1>())
    .def("__iter__", [](Table& self) {
        return nb::make_iterator(nb::type<Table>(), "iterator", self.begin(), self.end(),
                                 nb::keep_alive<0, 1>());
    }, nb::keep_alive<0, 1>())
    .def("__getitem__", &Table::at, nb::keep_alive<0, 1>())
    .def("__delitem__", &Table::remove_row)
    .def("__delitem__", &delitem_slice<Table>)
    .def("__bool__", &Table::ok)
    .def("__len__", &Table::length)
    .def("__repr__", [](const Table& self) {
        std::string s = "<gemmi.cif.Table ";
        if (self.ok())
          gemmi::cat_to(s, self.length(), " x ", self.width(), '>');
        else
          s += "nil>";
        return s;
    });

  cif_table_row
    .def_ro("row_index", &Table::Row::row_index)
    .def("__len__", &Table::Row::size)
    .def("__getitem__", (std::string& (Table::Row::*)(int)) &Table::Row::at)
    .def("__getitem__", [](Table::Row &self, const std::string& tag) {
        return self.value_at_unsafe(self.tab.find_column_position(tag));
    })
    .def("__setitem__", [](Table::Row &self, int idx, const std::string& value) {
        self.at(idx) = value;
    })
    .def("__setitem__", [](Table::Row &self, const std::string& tag,
                           const std::string& value) {
        self.value_at_unsafe(self.tab.find_column_position(tag)) = value;
    })
    .def("get", (std::string* (Table::Row::*)(int)) &Table::Row::ptr_at,
         nb::arg("index"), nb::rv_policy::reference_internal)
    .def("has", &Table::Row::has, nb::arg("index"))
    .def("__iter__", [](Table::Row& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>())
    .def("__repr__", [](const Table::Row& self) {
        std::string items;
        for (int i = 0; (size_t)i != self.size(); ++i)
          items += " " + (self.has(i) ? self[i] : "None");
        return "<gemmi.cif.Table.Row:" + items + ">";
    })
    .def("str", &Table::Row::str);

  cif.def("quote", &quote, nb::arg("string"));
  cif.def("quote_list", &quote_list);

  nb::class_<Ddl>(cif, "Ddl")
    .def("__init__", [](Ddl* ddl, gemmi::Logger&& logger,
                        bool print_unknown_tags, bool use_regex,
                        bool use_context, bool use_parents,
                        bool use_mandatory, bool use_unique_keys) {
          new(ddl) Ddl();
          ddl->logger = std::move(logger);
          ddl->print_unknown_tags = print_unknown_tags;
          ddl->use_regex = use_regex;
          ddl->use_context = use_context;
          ddl->use_parents = use_parents;
          ddl->use_mandatory = use_mandatory;
          ddl->use_unique_keys = use_unique_keys;
    }, nb::arg("logger"),
       nb::arg("print_unknown_tags")=true, nb::arg("use_regex")=true,
       nb::arg("use_context")=true, nb::arg("use_linked_groups")=true,
       nb::arg("use_mandatory")=true, nb::arg("use_unique_keys")=true)
    .def("set_logger", [](Ddl& self, gemmi::Logger&& logger) { self.logger = std::move(logger); })
    .def("read_ddl", [](Ddl& self, Document& doc) {
        self.read_ddl(std::move(doc));
        doc.clear();
    }, nb::arg("doc"))
    .def("validate_cif", &Ddl::validate_cif)
    ;
}
