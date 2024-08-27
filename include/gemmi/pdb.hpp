// Copyright 2017 Global Phasing Ltd.
//
// Read PDB file format and store it in Structure.
//
// Based on the format spec:
// https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
// + support for two-character chain IDs (columns 21 and 22)
// + read segment ID (columns 73-76)
// + read hybrid-36 serial numbers (http://cci.lbl.gov/hybrid_36/)
// + hybrid-36 sequence id for sequences longer than 9999 (no such examples)
// + allow for longer REMARK lines (up to 120 characters)

#ifndef GEMMI_PDB_HPP_
#define GEMMI_PDB_HPP_

#include <cstdio>     // for stdin, size_t
#include <unordered_map>
#include "fileutil.hpp" // for path_basename, file_open
#include "input.hpp"    // for FileStream
#include "model.hpp"    // for Structure, ...

namespace gemmi {

/// Returns operations corresponding to 1555, 2555, ... N555
GEMMI_DLL std::vector<Op> read_remark_290(const std::vector<std::string>& raw_remarks);

namespace impl {

struct GEMMI_DLL PdbReader {
  PdbReader(const PdbReadOptions& options_) : options(options_) {
    if (options.max_line_length <= 0 || options.max_line_length > 120)
      options.max_line_length = 120;
  }

  template<typename Stream>
  Structure from_stream(Stream&& stream, const std::string& source) {
    Structure st;
    st.input_format = CoorFormat::Pdb;
    st.name = path_basename(source, {".gz", ".pdb"});
    char line[122] = {0};
    while (size_t len = copy_line_from_stream(line, options.max_line_length+1, stream)) {
      ++line_num;
      read_pdb_line(line, len, st, source);
      if (is_end)
        break;
    }
    finalize_structure_after_reading_pdb(st);
    return st;
  }

private:
  int line_num = 0;
  bool after_ter = false;
  bool is_end = false;
  PdbReadOptions options;
  Model *model = nullptr;
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  Transform matrix;
  std::vector<std::string> conn_records;
  std::unordered_map<ResidueId, int> resmap;

  [[noreturn]] void wrong(const std::string& msg) const {
    fail("Problem in line ", std::to_string(line_num), ": " + msg);
  }
  void read_pdb_line(const char* line, size_t len, Structure& st, const std::string& source);
  void finalize_structure_after_reading_pdb(Structure& st) const;
};

}  // namespace impl

inline Structure read_pdb_file(const std::string& path,
                               PdbReadOptions options=PdbReadOptions()) {
  auto f = file_open(path.c_str(), "rb");
  return impl::PdbReader(options).from_stream(FileStream{f.get()}, path);
}

inline Structure read_pdb_from_memory(const char* data, size_t size,
                                      const std::string& name,
                                      PdbReadOptions options=PdbReadOptions()) {
  return impl::PdbReader(options).from_stream(MemoryStream(data, size), name);
}

inline Structure read_pdb_string(const std::string& str,
                                 const std::string& name,
                                 PdbReadOptions options=PdbReadOptions()) {
  return read_pdb_from_memory(str.c_str(), str.length(), name, options);
}

// A function for transparent reading of stdin and/or gzipped files.
template<typename T>
inline Structure read_pdb(T&& input, PdbReadOptions options=PdbReadOptions()) {
  if (input.is_stdin())
    return impl::PdbReader(options).from_stream(FileStream{stdin}, "stdin");
  if (input.is_compressed())
    return impl::PdbReader(options).from_stream(input.get_uncompressing_stream(), input.path());
  return read_pdb_file(input.path(), options);
}

} // namespace gemmi
#endif
