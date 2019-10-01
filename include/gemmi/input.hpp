// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.
// Used to decouple file reading and uncompression.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <cassert>
#include <cstring> // for memchr
#include <memory>  // for unique_ptr
#include <string>
#include "fail.hpp"  // for unreachable

namespace gemmi {

// provides the same methods as GzStream
struct FileStream {
  std::FILE* f;
  // used in pdb.hpp
  char* gets(char* line, int size) { return std::fgets(line, size, f); }
  int getc() { return std::fgetc(f); }
  // used in ccp4.hpp
  bool read(void* buf, size_t len) { return std::fread(buf, len, 1, f) == 1; }
  bool seek(long offset) { return std::fseek(f, offset, SEEK_SET) == 0; }
};

struct MemoryStream {
  MemoryStream(const char* start_, const char* end_)
    : start(start_), end(end_), cur(start_) {}

  char* gets(char* line, int size) {
    if (cur >= end)
      return nullptr;
    if (size > end - cur)
      size = end - cur;
    const char* nl = (const char*) std::memchr(cur, '\n', size);
    size_t len = nl ? nl - cur + 1 : size;
    std::memcpy(line, cur, len);
    cur += len;
    return line;
  }
  int getc() { return cur < end ? *++cur : EOF; }

  bool read(void* buf, size_t len) {
    if (cur + len > end)
      return false;
    std::memcpy(buf, cur, len);
    cur += len;
    return true;
  }

  int seek(long offset) {
    cur = start + offset;
    return cur < end;
  }

private:
  const char* const start;
  const char* const end;
  const char* cur;
};


class BasicInput {
public:
  explicit BasicInput(const std::string& path) : path_(path) {}

  const std::string& path() const { return path_; };
  const std::string& basepath() const { return path_; };

  // Does the path stands for stdin?
  // Each reading function needs to call it (some functions use stdin
  // and some std::cin, so we don't try to unify it here).
  bool is_stdin() const { return path() == "-"; };

  // providing the same interface as MaybeGzipped
  bool is_compressed() const { return false; }
  FileStream get_uncompressing_stream() const { assert(0); unreachable(); }
  // for reading (uncompressing into memory) the whole file at once
  std::unique_ptr<char[]> memory() { return nullptr; }
  size_t memory_size() const { return 0; };

private:
  std::string path_;
};

} // namespace gemmi
#endif
