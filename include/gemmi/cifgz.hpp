// Copyright 2017 Global Phasing Ltd.
//
// Functions for transparent reading of gzipped files. Uses zlib.

#ifndef GEMMI_CIFGZ_HPP_
#define GEMMI_CIFGZ_HPP_
#include "cif.hpp"
#include "util.hpp"  // ends_with
#include <cstdio>
#include <zlib.h>

namespace gemmi {
namespace cif {

// Throws if the size is not found or if it is suspicious.
// Anything outside of the arbitrary limits from 1 to 10x of the compressed
// size looks suspicious to us.
inline size_t estimate_uncompressed_size(const std::string& path) {
  std::FILE* file;
#ifdef _WIN32
  if (::fopen_s(&file, path.c_str(), "rb") != 0)
#else
  if ((file = std::fopen(path.c_str(), "rb")) == nullptr)
#endif
    fail("Failed to open file: " + path);
  std::unique_ptr<FILE, decltype(&std::fclose)> cleanup(file, &std::fclose);
  if (std::fseek(file, -4, SEEK_END) != 0)
    fail("fseek() failed (empty file?): " + path);
  long pos = std::ftell(file);
  if (pos <= 0)
    fail("ftell() failed on " + path);
  size_t gzipped_size = pos + 4;
  unsigned char buf[4];
  if (std::fread(buf, 1, 4, file) != 4)
    fail("Failed to read last 4 bytes of: " + path);
  size_t orig_size = (buf[3] << 24) | (buf[2] << 16) | (buf[1] << 8) | buf[0];
  if (orig_size < gzipped_size || orig_size > 10 * gzipped_size)
    fail("Cannot determine uncompressed size of " + path);
  return orig_size;
}


inline std::unique_ptr<char[]> gunzip_to_memory(const std::string& path,
                                                size_t orig_size) {
  if (orig_size > 500000000)
    fail("For now gz files above 500MB uncompressed are not supported.");
  std::unique_ptr<char[]> mem(new char[orig_size]);
  gzFile file = gzopen(path.c_str(), "rb");
  if (!file)
    fail("Failed to gzopen: " + path);
  int bytes_read = gzread(file, mem.get(), orig_size);
  if (bytes_read < (int) orig_size && !gzeof(file)) {
    int errnum;
    std::string err_str = gzerror(file, &errnum);
    if (errnum) {
      gzclose(file);
      fail("Error reading " + path + ": " + err_str);
    }
  }
  gzclose(file);
  return mem;
}

inline Document read_any(const std::string& path) {
  if (path == "-")
    return read_cstream(stdin, 16*1024, "stdin");
  if (ends_with(path, ".gz")) {
    size_t orig_size = estimate_uncompressed_size(path);
    std::unique_ptr<char[]> mem = cif::gunzip_to_memory(path, orig_size);
    return read_memory(mem.get(), orig_size, path.c_str());
  }
  return read_file(path);
}

inline bool check_syntax_any(const std::string& path, std::string* msg) {
  if (gemmi::ends_with(path, ".gz")) {
    size_t orig_size = cif::estimate_uncompressed_size(path);
    std::unique_ptr<char[]> mem = cif::gunzip_to_memory(path, orig_size);
    tao::pegtl::memory_input<> in(mem.get(), orig_size, path);
    return cif::check_syntax(in, msg);
  }
  tao::pegtl::file_input<> in(path);
  return cif::check_syntax(in, msg);
}

} // namespace cif
} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
