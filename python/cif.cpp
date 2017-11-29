// Copyright 2017 Global Phasing Ltd.

#include "gemmi/cifdoc.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi::cif;

void add_cif(py::module& cif) {
  py::class_<Document>(cif, "Document")
    .def(py::init<>())
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
    .def("__delitem__", [](Document &d, int index) {
        if (index < 0)
          index += d.blocks.size();
        if (index < 0 || static_cast<size_t>(index) >= d.blocks.size())
          throw py::index_error();
        d.blocks.erase(d.blocks.begin() + index);
    }, py::arg("index"))
    .def("clear", &Document::clear)
    .def("sole_block", &Document::sole_block,
         py::return_value_policy::reference_internal,
         "Returns the only block if there is exactly one")
    .def("find_block", &Document::find_block, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("write_file", &write_to_file, py::arg("filename"),
         "Write data to a CIF file.")
    .def("as_string", [](const Document& d) {
        std::ostringstream os;
        os << d;
        return os.str();
    }, "Write data in a CIF format to a string.")
    .def("as_json", [](const Document& d) {
        std::ostringstream os;
        JsonWriter(os).write_json(d);
        return os.str();
    }, "Returns JSON representation in a string.");

  py::class_<Block>(cif, "Block")
    .def(py::init<>())
    .def_readwrite("name", &Block::name)
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
         py::arg("prefix"), py::arg("tags"))
    .def("find", (Table (Block::*)(const std::vector<std::string>&))
                 &Block::find,
         py::arg("tags"))
    .def("set_pair", &Block::set_pair, py::arg("tag"), py::arg("value"))
    .def("init_loop", &Block::init_loop, py::arg("prefix"), py::arg("tags"),
         py::return_value_policy::reference_internal)
    .def("find_mmcif_category", &Block::find_mmcif_category,
         py::arg("category"), py::keep_alive<0, 1>(),
         "Returns Table with all items in the category.")
    .def("get_mmcif_category_names", &Block::get_mmcif_category_names,
         "For mmCIF files only. Returns list of all category prefixes (_x.)")
    .def("init_mmcif_loop", &Block::init_mmcif_loop,
         py::arg("cat"), py::arg("tags"),
         py::return_value_policy::reference_internal)
    .def("__repr__", [](const Block &self) {
        return "<gemmi.cif.Block " + self.name + ">";
    });

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
        return "<gemmi.cif.Loop " + std::to_string(self.length()) + " x " +
                                    std::to_string(self.width()) + ">";
    });

  py::class_<Column>(cif, "Column")
    .def(py::init<>())
    .def("get_loop", &Column::get_loop,
         py::return_value_policy::reference_internal)
    .def("__iter__", [](const Column& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    .def("__bool__", [](const Column &self) -> bool { return self.item(); })
    .def("__getitem__", (std::string& (Column::*)(int)) &Column::at)
    .def("__setitem__", [](Column &self, int idx, std::string value) {
        self.at(idx) = value;
    })
    .def("str", &Column::str)
    .def("__repr__", [](const Column &self) {
        std::string desc = "nil";
        if (const std::string* tag = self.get_tag())
          desc = *tag + " length " + std::to_string(self.length());
        return "<gemmi.cif.Column " + desc + ">";
    });

  py::class_<Table> lt(cif, "Table");
  lt.def("width", &Table::width)
    .def("column", &Table::column, py::arg("n"), py::keep_alive<0, 1>())
    .def("find_row", &Table::find_row, py::keep_alive<0, 1>())
    .def("find_column", &Table::find_column, py::arg("suffix"),
         py::keep_alive<0, 1>())
    .def("erase", &Table::erase)
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
               (self.ok() ? std::to_string(self.length()) + " x " +
                            std::to_string(self.width())
                          : "nil") +
               ">";
    });

  py::class_<Table::Row>(lt, "Row")
    .def("str", &Table::Row::str)
    .def("__len__", &Table::Row::size)
    .def("__getitem__", (std::string& (Table::Row::*)(int)) &Table::Row::at)
    .def("__setitem__", [](Table::Row &self, int idx, std::string value) {
        self.at(idx) = value;
    })
    .def("get", &Table::Row::ptr_at,
         py::arg("index"), py::return_value_policy::reference_internal)
    .def("__iter__", [](const Table::Row& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [](const Table::Row& self) {
        std::string items;
        for (size_t i = 0; i != self.size(); ++i)
          items += " " + (self.has(i) ? self[i] : "None");
        return "<gemmi.cif.Table.Row:" + items + ">";
    });
}
