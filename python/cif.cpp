// Copyright 2017 Global Phasing Ltd.

#include "gemmi/cifdoc.hpp"
#include "gemmi/tostr.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/ofstream.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi::cif;

void add_cif(py::module& cif) {
  py::enum_<Style>(cif, "Style")
    .value("Simple", Style::Simple)
    .value("NoBlankLines", Style::NoBlankLines)
    .value("PreferPairs", Style::PreferPairs)
    .value("Pdbx", Style::Pdbx)
    .value("Indent35", Style::Indent35);
  py::class_<Document>(cif, "Document")
    .def(py::init<>())
    .def_readwrite("source", &Document::source)
    .def("__len__", [](const Document& d) { return d.blocks.size(); })
    .def("__iter__", [](const Document& d) {
        return py::make_iterator(d.blocks);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Document& d, const std::string& name) -> Block& {
        Block* b = d.find_block(name);
        if (!b)
          throw py::key_error("block '" + name + "' does not exist");
        return *b;
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("__getitem__", [](Document& d, int index) -> Block& {
        return d.blocks.at(index >= 0 ? index : index + d.blocks.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", [](Document &d, py::slice slice) -> py::list {
        size_t start, stop, step, slicelength;
        if (!slice.compute(d.blocks.size(), &start, &stop, &step, &slicelength))
          throw py::error_already_set();
        py::list l;
        for (size_t i = 0; i < slicelength; ++i)
          l.append(py::cast(d.blocks[start + i * step]));
        return l;
    })
    .def("__delitem__", [](Document &d, int index) {
        if (index < 0)
          index += (int) d.blocks.size();
        if (index < 0 || static_cast<size_t>(index) >= d.blocks.size())
          throw py::index_error();
        d.blocks.erase(d.blocks.begin() + index);
    }, py::arg("index"))
    .def("add_copied_block", [](Document& d, const Block& block, int pos) {
        if (pos < 0)
          pos = (int) d.blocks.size();
        else if (pos > (int) d.blocks.size())
          throw py::index_error();
        d.blocks.insert(d.blocks.begin() + pos, block);
        return d.blocks[pos];
    }, py::arg("block"), py::arg("pos")=-1,
       py::return_value_policy::reference_internal)
    .def("add_new_block", &Document::add_new_block,
         py::arg("name"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("clear", &Document::clear)
    .def("sole_block", (Block& (Document::*)()) &Document::sole_block,
         py::return_value_policy::reference_internal,
         "Returns the only block if there is exactly one")
    .def("find_block",
         (Block* (Document::*)(const std::string&)) &Document::find_block,
         py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("write_file",
         [](const Document& doc, const std::string& filename, Style s) {
        gemmi::Ofstream os(filename);
        write_cif_to_stream(os.ref(), doc, s);
    }, py::arg("filename"), py::arg("style")=Style::Simple,
    "Write data to a CIF file.")
    .def("as_string", [](const Document& d, Style style) {
        std::ostringstream os;
        write_cif_to_stream(os, d, style);
        return os.str();
    }, py::arg("style")=Style::Simple, "Returns a string in CIF format.")
    .def("as_json", [](const Document& d, bool mmjson) {
        std::ostringstream os;
        JsonWriter writer(os);
        if (mmjson)
          writer.set_mmjson();
        writer.write_json(d);
        return os.str();
    }, py::arg("mmjson")=false, "Returns JSON representation in a string.")
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

  py::class_<Block>(cif, "Block")
    .def(py::init<const std::string &>())
    .def_readwrite("name", &Block::name)
    .def("__iter__", [](Block& self) { return py::make_iterator(self.items); },
         py::keep_alive<0, 1>())
    .def("find_pair", &Block::find_pair, py::arg("tag"),
         py::return_value_policy::reference_internal)
    .def("find_value", &Block::find_value, py::arg("tag"),
         py::return_value_policy::reference)
    .def("find_loop", &Block::find_loop, py::arg("tag"),
         py::keep_alive<0, 1>())
    .def("find_values", &Block::find_values, py::arg("tag"),
         py::keep_alive<0, 1>())
    .def("find", (Table (Block::*)(const std::string&,
            const std::vector<std::string>&)) &Block::find,
         py::arg("prefix"), py::arg("tags"), py::keep_alive<0, 1>())
    .def("find", (Table (Block::*)(const std::vector<std::string>&))
                 &Block::find,
         py::arg("tags"), py::keep_alive<0, 1>())
    .def("find_or_add", &Block::find_or_add,
         py::arg("prefix"), py::arg("tags"), py::keep_alive<0, 1>())
    .def("find_frame", &Block::find_frame, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("get_index", &Block::get_index, py::arg("tag"))
    .def("set_pair", &Block::set_pair, py::arg("tag"), py::arg("value"))
    .def("init_loop", &Block::init_loop, py::arg("prefix"), py::arg("tags"),
         py::return_value_policy::reference_internal)
    .def("move_item", &Block::move_item, py::arg("old_pos"), py::arg("new_pos"))
    .def("find_mmcif_category", &Block::find_mmcif_category,
         py::arg("category"), py::keep_alive<0, 1>(),
         "Returns Table with all items in the category.")
    .def("get_mmcif_category_names", &Block::get_mmcif_category_names,
         "For mmCIF files only. Returns list of all category prefixes (_x.)")
    .def("init_mmcif_loop", &Block::init_mmcif_loop,
         py::arg("cat"), py::arg("tags"),
         py::return_value_policy::reference_internal)
    .def("set_mmcif_category",
         [](Block &self, std::string name, py::dict data, bool raw) {
           size_t w = data.size();
           std::vector<std::string> tags;
           tags.reserve(w);
           std::vector<py::list> values;
           values.reserve(w);
           for (auto item : data) {
             tags.emplace_back(py::str(item.first));
             values.emplace_back(item.second.cast<py::list>());
             if (values.back().size() != values[0].size())
               throw py::value_error("all columns must have equal length");
           }
           if (w == 0 || values[0].size() == 0)
             throw py::value_error("data cannot be empty");
           Loop& loop = self.init_mmcif_loop(std::move(name), std::move(tags));
           loop.values.resize(w * values[0].size());
           for (size_t col = 0; col != w; ++col) {
             size_t idx = col;
             for (auto handle : values[col]) {
               std::string& val = loop.values[idx];
               PyObject* ptr = handle.ptr();
               if (handle.is_none()) {
                 val = "?";
               } else if (ptr == Py_False) {
                 val = ".";
               } else if (ptr == Py_True) {
                 throw py::value_error("unexpected value True");
               } else if (raw || PyFloat_Check(ptr) || PyLong_Check(ptr)) {
                 val = py::str(handle);
               } else {
                 val = quote(py::str(handle));
               }
               idx += w;
             }
           }
         }, py::arg("name"), py::arg("data"), py::arg("raw")=false)
    .def("get_mmcif_category",
         [](Block &self, std::string name, bool raw) {
           ensure_mmcif_category(name);
           py::dict data;
           Table table = self.find_mmcif_category(name);
           int len = (int) table.length();
           for (const std::string& tag : table.tags()) {
             assert(tag.size() >= name.size());
             Column col = table.column((int)data.size());
             py::list new_list(len);
             for (int i = 0; i != len; ++i)
               if (raw) {
                 new_list[i] = col[i];
               } else if (is_null(col[i])) {
                 if (col[i][0] == '?')
                   new_list[i] = py::none().inc_ref();
                 else
                   new_list[i] = py::handle(Py_False).inc_ref();
               } else {
                 new_list[i] = col.str(i);
               }
             data[tag.c_str() + name.size()] = new_list;
           }
           return data;
         }, py::arg("name"), py::arg("raw")=false)
    .def("__repr__", [](const Block &self) {
        return gemmi::tostr("<gemmi.cif.Block ", self.name, '>');
    });

  py::class_<Item> (cif, "Item")
    .def("erase", &Item::erase)
    .def_readonly("line_number", &Item::line_number)
    .def_property_readonly("pair", [](Item& self) {
        return self.type == ItemType::Pair ? &self.pair : nullptr;
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("loop", [](Item& self) {
        return self.type == ItemType::Loop ? &self.loop : nullptr;
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("frame", [](Item& self) {
        return self.type == ItemType::Frame ? &self.frame : nullptr;
    }, py::return_value_policy::reference_internal)
    ;

  py::class_<Loop> lp(cif, "Loop");
  lp.def(py::init<>())
    .def("width", &Loop::width, "Returns number of columns")
    .def("length", &Loop::length, "Returns number of rows")
    .def_readonly("tags", &Loop::tags)
    .def("val", &Loop::val, py::arg("row"), py::arg("col"))
    .def("add_row", &Loop::add_row<std::vector<std::string>>,
         py::arg("new_values"), py::arg("pos")=-1)
    .def("set_all_values", &Loop::set_all_values, py::arg("columns"))
    .def("__repr__", [](const Loop &self) {
        return gemmi::tostr("<gemmi.cif.Loop ", self.length(), " x ",
                            self.width(), '>');
    });

  py::class_<Column>(cif, "Column")
    .def(py::init<>())
    .def("get_loop", &Column::get_loop,
         py::return_value_policy::reference_internal)
    .def("__iter__", [](const Column& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    .def("__bool__", [](const Column &self) { return self.item() != nullptr; })
    .def("__len__", [](const Column &self) { return self.length(); })
    .def("__getitem__", (std::string& (Column::*)(int)) &Column::at)
    .def("__setitem__", [](Column &self, int idx, std::string value) {
        self.at(idx) = value;
    })
    .def("str", &Column::str, py::arg("index"))
    .def("__repr__", [](const Column &self) {
        std::string desc = "nil";
        if (const std::string* tag = self.get_tag())
          desc = gemmi::tostr(*tag, " length ", self.length());
        return "<gemmi.cif.Column " + desc + ">";
    });

  py::class_<Table> lt(cif, "Table");
  lt.def("width", &Table::width)
    .def_readonly("prefix_length", &Table::prefix_length)
    .def_property_readonly("loop", &Table::get_loop,
                           py::return_value_policy::reference_internal)
    .def("get_prefix", &Table::get_prefix)
    .def("has_column", &Table::has_column)
    .def("column", &Table::column, py::arg("n"), py::keep_alive<0, 1>())
    .def("find_row", &Table::find_row, py::keep_alive<0, 1>())
    .def("find_column", &Table::find_column, py::arg("tag"),
         py::keep_alive<0, 1>())
    .def("erase", &Table::erase)
    .def("append_row", &Table::append_row<std::vector<std::string>>,
         py::arg("new_values"))
    .def_property_readonly("tags",
            py::cpp_function(&Table::tags, py::keep_alive<0, 1>()))
    .def("__iter__", [](Table& self) {
        return py::make_iterator(self, py::keep_alive<0, 1>());
    }, py::keep_alive<0, 1>())
    .def("__getitem__", &Table::at, py::keep_alive<0, 1>())
    .def("__bool__", &Table::ok)
    .def("__len__", &Table::length)
    .def("__repr__", [](const Table& self) {
        return "<gemmi.cif.Table " +
               (self.ok() ? gemmi::tostr(self.length(), " x ", self.width())
                          : "nil") +
               ">";
    });

  py::class_<Table::Row>(lt, "Row")
    .def("str", &Table::Row::str)
    .def("__len__", &Table::Row::size)
    .def("__getitem__", (std::string& (Table::Row::*)(int)) &Table::Row::at)
    .def("__getitem__", [](Table::Row &self, const std::string& tag) {
        return self.value_at_unsafe(self.tab.find_column_position(tag));
    })
    .def("__setitem__", [](Table::Row &self, int idx, std::string value) {
        self.at(idx) = value;
    })
    .def("__setitem__", [](Table::Row &self, const std::string& tag,
                           std::string value) {
        self.value_at_unsafe(self.tab.find_column_position(tag)) = value;
    })
    .def("get", (std::string* (Table::Row::*)(int)) &Table::Row::ptr_at,
         py::arg("index"), py::return_value_policy::reference_internal)
    .def("has", &Table::Row::has, py::arg("index"))
    .def("__iter__", [](const Table::Row& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [](const Table::Row& self) {
        std::string items;
        for (size_t i = 0; i != self.size(); ++i)
          items += " " + (self.has(i) ? self[i] : "None");
        return "<gemmi.cif.Table.Row:" + items + ">";
    });

  cif.def("quote", &quote, py::arg("string"));
}
