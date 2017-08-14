// Copyright 2017 Global Phasing Ltd.
//
// Functions for transparent reading of gzipped files. Uses zlib.

#ifndef GEMMI_GZ_HPP_
#define GEMMI_GZ_HPP_
#include "util.hpp"  // ends_with, MaybeStdin
#include <cstdio>
#include <zlib.h>

namespace gemmi {

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
  std::unique_ptr<FILE, decltype(&fclose)> cleanup(file, &fclose);
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

class MaybeGzipped : public MaybeStdin {
public:
  explicit MaybeGzipped(const std::string& path)
    : MaybeStdin(path), mem_size_(0) {}
  size_t mem_size() const { return mem_size_; };
  std::unique_ptr<char[]> memory() {
    if (!ends_with(path(), ".gz"))
      return MaybeStdin::memory();
    mem_size_ = estimate_uncompressed_size(path());
    return gunzip_to_memory(path(), mem_size_);
  }
private:
  size_t mem_size_;
};

} // namespace gemmi

#endif
// vim:sw=2:ts=2:et
