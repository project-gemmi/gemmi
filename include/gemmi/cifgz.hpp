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
    throw std::runtime_error("Failed to open file: " + path);
  std::unique_ptr<FILE, decltype(&fclose)> cleanup(file, &fclose);
  if (std::fseek(file, -4, SEEK_END) != 0)
    throw std::runtime_error("fseek() failed (empty file?): " + path);
  long pos = std::ftell(file);
  if (pos <= 0)
    throw std::runtime_error("ftell() failed on " + path);
  size_t gzipped_size = pos + 4;
  unsigned char buf[4];
  if (std::fread(buf, 1, 4, file) != 4)
    throw std::runtime_error("Failed to read last 4 bytes of: " + path);
  size_t orig_size = (buf[3] << 24) | (buf[2] << 16) | (buf[1] << 8) | buf[0];
  if (orig_size < gzipped_size || orig_size > 10 * gzipped_size)
    throw std::runtime_error("Cannot determine uncompressed size of " + path);
  return orig_size;
}


inline std::unique_ptr<char[]> gunzip_to_memory(const std::string& path,
                                                size_t orig_size) {
  if (orig_size > 500000000)
    throw std::runtime_error("For now gz files above 500MB uncompressed"
                             " are not supported.");
  std::unique_ptr<char[]> mem(new char[orig_size]);
  gzFile file = gzopen(path.c_str(), "rb");
  if (!file)
    throw std::runtime_error("Failed to gzopen: " + path);
  int bytes_read = gzread(file, mem.get(), orig_size);
  if (bytes_read < (int) orig_size && !gzeof(file)) {
    int errnum;
    std::string err_str = gzerror(file, &errnum);
    if (errnum) {
      gzclose(file);
      throw std::runtime_error("Error reading " + path + ": " + err_str);
    }
  }
  gzclose(file);
  return mem;
}

inline Document read_any(const std::string& path) {
  Document d;
  if (path == "-") {
    d.read_cstream(stdin, 16*1024, "stdin");
  } else if (ends_with(path, ".gz")) {
    size_t orig_size = estimate_uncompressed_size(path);
    std::unique_ptr<char[]> mem = cif::gunzip_to_memory(path, orig_size);
    d.read_memory(mem.get(), orig_size, path.c_str());
  } else {
    d.read_file(path);
  }
  return d;
};

} // namespace cif
} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
