// Copyright 2017 Global Phasing Ltd.
//
// Functions for transparent reading of gzipped files. Uses zlib.

#ifndef GEMMI_CIFGZ_HH_
#define GEMMI_CIFGZ_HH_
#include "cif.hh"
#include <stdio.h>
#include <zlib.h>

namespace gemmi {
namespace cif {

// Returns 0 if the size is not found or if it is suspicious.
// Anything outside of the arbitrary limits from 1 to 10x of the compressed
// size looks suspicious to us.
inline size_t estimate_uncompressed_size(const std::string& path) {
  std::FILE* file;
#ifdef _WIN32
  if (::fopen_s(&file, path.c_str(), "rb") != 0)
#else
  if ((file = std::fopen(path.c_str(), "rb")) == nullptr)
#endif
    return 0;
  std::unique_ptr<FILE, decltype(&fclose)> cleanup(file, &fclose);
  if (std::fseek(file, -4, SEEK_END) != 0)
    return 0; // fseek() failed
  long pos = std::ftell(file);
  if (pos <= 0)
    return 0; // ftell() failed
  size_t gzipped_size = pos + 4;
  unsigned char buf[4];
  if (std::fread(buf, 1, 4, file) != 4)
    return 0; // fread() failed
  size_t orig_size = (buf[3] << 24) | (buf[2] << 16) | (buf[1] << 8) | buf[0];
  if (orig_size < gzipped_size || orig_size > 10 * gzipped_size)
    return 0;
  return orig_size;
}


inline void gunzip_and_read(Document& d, const std::string& path) {
  size_t orig_size = estimate_uncompressed_size(path);
  if (orig_size == 0)
    throw std::runtime_error("Failed to open or read file: " + path);
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
  d.read_memory(mem.get(), orig_size, path.c_str());
}

inline Document read_any(const std::string& path) {
  Document d;
  if (path == "stdin") { // temporary, Clara can't handle "-"
    d.read_cstream(stdin, "stdin", 16*1024);
    //d.read_istream(std::cin, "stdin", 16*1024);
  } else if (path.size() > 3 && path.substr(path.size() - 3) == ".gz") {
    cif::gunzip_and_read(d, path);
  } else {
    d.read_file(path);
  }
  return d;
};

} // namespace cif
} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
