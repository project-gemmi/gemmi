// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.
// Used to decouple file reading and uncompression.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <cassert>
#include <cstring> // for memchr
#include <cstdlib> // for malloc, realloc
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
  MemoryStream(const char* start_, size_t size)
    : start(start_), end(start_ + size), cur(start_) {}

  char* gets(char* line, int size) {
    if (cur >= end)
      return nullptr;
    if (size > end - cur)
      size = int(end - cur);
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

class CharArray {
  std::unique_ptr<char[], decltype(std::free)*> ptr_;
  size_t size_;
public:
  CharArray() : ptr_(nullptr, std::free), size_(0) {}
  explicit CharArray(size_t n) : ptr_((char*)std::malloc(n), std::free), size_(n) {};
  explicit operator bool() const { return (bool)ptr_; }
  char* data() { return ptr_.get(); }
  const char* data() const { return ptr_.get(); }
  size_t size() const { return size_; }
  void set_size(size_t n) { size_ = n; }

  void resize(size_t n) {
    char* new_ptr = (char*) std::realloc(ptr_.get(), n);
    if (!new_ptr)
      fail("Out of memory.");
    ptr_.release();
    ptr_.reset(new_ptr);
    size_ = n;
  }
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
  CharArray uncompress_into_buffer() { return {}; }

private:
  std::string path_;
};

} // namespace gemmi
#endif
