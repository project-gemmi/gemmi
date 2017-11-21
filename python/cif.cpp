// Copyright 2017 Global Phasing Ltd.

#include "gemmi/numb.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;
using namespace gemmi::cif;

void init_cif(py::module& cif) {
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
         "Returns the only block if there is exactly one")
    .def("find_block", &Document::find_block, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("write_file", &write_to_file, py::arg("filename"),
         "Write data to a CIF file.")
    .def("as_json", [](const Document& d) {
        std::ostringstream os;
        JsonWriter(os).write_json(d);
        return os.str();
    }, "Returns JSON representation in a string.");

  py::class_<Block>(cif, "Block")
    .def(py::init<>())
    .def_readonly("name", &Block::name)
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
    .def("find_mmcif_category", &Block::find_mmcif_category,
         py::arg("category"), py::keep_alive<0, 1>(),
         "Returns Table with all items in the category.")
    .def("delete_category", &Block::delete_category, py::arg("prefix"),
         "End mmCIF category with the dot: block.delete_category('_exptl.')")
    .def("get_mmcif_category_names", &Block::get_mmcif_category_names,
         "For mmCIF files only. Returns list of all category prefixes (_x.)")
    .def("__repr__", [](const Block &self) {
        return "<gemmi.cif.Block " + self.name + ">";
    });

  py::class_<Loop> lp(cif, "Loop");
  lp.def(py::init<>())
    .def("width", &Loop::width, "Returns number of columns")
    .def("length", &Loop::length, "Returns number of rows")
    .def_readonly("tags", &Loop::tags)
    .def("val", &Loop::val, py::arg("row"), py::arg("col"))
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
    .def("str", &Column::str)
    .def("__repr__", [](const Column &self) {
        std::string desc = "nil";
        if (const std::string* tag = self.get_tag())
          desc = *tag + " length " + std::to_string(self.length());
        return "<gemmi.cif.Column " + desc + ">";
    });

  py::class_<Table> lt(cif, "Table");
  lt.def_readonly("loop", &Table::loop)
    .def("find_row", &Table::find_row, py::keep_alive<0, 1>())
    .def("__iter__", [](Table& self) {
        return py::make_iterator(self, py::keep_alive<0, 1>());
    }, py::keep_alive<0, 1>())
    .def("__getitem__", &Table::at, py::keep_alive<0, 1>())
    .def("__bool__", &Table::ok)
    .def("__len__", &Table::length)
    .def("width", &Table::width)
    .def_property_readonly("tags",
            py::cpp_function(&Table::tags, py::keep_alive<0, 1>()))
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

  cif.def("read_file", &read_file, py::arg("filename"),
          "Reads a CIF file copying data into Document.");
  cif.def("read", [](const std::string& s) { return read(MaybeGzipped(s)); },
          py::arg("filename"), "Reads normal or gzipped CIF file.");
  cif.def("read_string", &read_string, py::arg("data"),
          "Reads a string as a CIF file.");

  cif.def("as_string", (std::string (*)(const std::string&)) &as_string,
          py::arg("value"), "Get string content (no quotes) from raw string.");
  cif.def("as_number", &as_number, py::arg("value"), py::arg("default")=NAN,
          "Returns float number from string");
  cif.def("as_int", (int (*)(const std::string&)) &as_int, py::arg("value"),
          "Returns int number from string value.");
  cif.def("as_int", (int (*)(const std::string&, int)) &as_int,
          py::arg("value"), py::arg("default"),
          "Returns int number from string value or the second arg if null.");
}
