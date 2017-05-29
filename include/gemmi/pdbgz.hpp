// Copyright 2017 Global Phasing Ltd.
//
// Functions for transparent reading of gzipped PDB files. Uses zlib.

#ifndef GEMMI_PDBGZ_HPP_
#define GEMMI_PDBGZ_HPP_
#include "pdb.hpp"
#include "util.hpp"  // ends_with
#include <cstdio>
#include <zlib.h>

namespace gemmi {
namespace mol {

namespace internal {

class GzipLineInput {
public:
  std::string source;
  GzipLineInput(std::string path) : source(path) {
    f_ = gzopen(source.c_str(), "rb");
    if (!f_)
      throw std::runtime_error("Failed to open file: " + path);
    gzbuffer(f_, 64*1024);
  }
  ~GzipLineInput() { gzclose(f_); }

  size_t copy_line(char* line) {
    if (!gzgets(f_, line, 82))
      return 0;
    size_t len = std::strlen(line);
    // We don't expect lines longer than 80 characters, but if one is found,
    // just discard the rest of the line.
    if (len > 0 && line[len-1] != '\n')
      for (int c = gzgetc(f_); c != 0 && c != -1 && c != '\n'; c = gzgetc(f_))
        continue;
    ++line_num_;
    return len;
  }
  size_t line_num() const { return line_num_; }
private:
  gzFile f_;
  size_t line_num_ = 0;
};

} // namespace internal

inline Structure read_pdb_gz(const std::string& path) {
  internal::GzipLineInput input(path);
  return read_pdb_from_input(input);
}

inline Structure read_pdb_any(const std::string& path) {
  if (path == "-")
    return read_pdb_from_cstream(stdin, "stdin");
  if (ends_with(path, ".gz"))
    return read_pdb_gz(path);
  return read_pdb(path);
}

} // namespace mol
} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
