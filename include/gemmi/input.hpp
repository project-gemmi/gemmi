// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.
// Used to decouple file reading and decompression.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <cassert>
#include <cstddef> // for ptrdiff_t
#include <cstdio>  // for FILE, fseek, fread
#include <cstdlib> // for malloc, realloc
#include <cstring> // for memchr
#include <string>
#include "fileutil.hpp"  // for fileptr_t

namespace gemmi {

// base class for FileStream, MemoryStream and GzStream
struct AnyStream {
  virtual ~AnyStream() = default;

  virtual char* gets(char* line, int size) = 0;   // for pdb, copy_line()
  virtual int getc() = 0;                         // for copy_line()
  virtual bool read(void* buf, size_t len) = 0;   // for ccp4, mtz

  // these are not used in GzStream because MemoryStream is used for mtz
  virtual bool seek(std::ptrdiff_t /*offset*/) { return false; } // for mtz
  virtual std::string read_rest() { return {}; }  // for mtz

  size_t copy_line(char* line, int size) {        // for pdb, xds_ascii
    if (!gets(line, size))
      return 0;
    size_t len = std::strlen(line);
    // If a line is longer than size we discard the rest of it.
    if (len > 0 && line[len-1] != '\n')
      for (int c = getc(); c > 0 /* not 0 nor EOF */ && c != '\n'; c = getc())
        continue;
    return len;
  };
};

struct FileStream final : public AnyStream {
  FileStream(std::FILE* f_) : f(f_, needs_fclose{false}) {}
  FileStream(const char* path, const char* mode) : f(file_open_or(path, mode, stdin)) {}

  char* gets(char* line, int size) override { return std::fgets(line, size, f.get()); }
  int getc() override { return std::fgetc(f.get()); }
  bool read(void* buf, size_t len) override { return std::fread(buf, len, 1, f.get()) == 1; }

  std::string read_rest() override {
    std::string ret;
    int c = std::fgetc(f.get());
    if (c != EOF) {
      ret += (char)c;
      char buf[512];
      for (;;) {
        size_t n = std::fread(buf, 1, sizeof(buf), f.get());
        ret.append(buf, n);
        if (n != sizeof(buf))
          break;
      }
    }
    return ret;
  }

  bool seek(std::ptrdiff_t offset) override {
#if defined(_MSC_VER)
    return _fseeki64(f.get(), offset, SEEK_SET) == 0;
#elif defined(__MINGW32__)
    return fseeko(f.get(), (_off_t)offset, SEEK_SET) == 0;
#else
    return std::fseek(f.get(), (long)offset, SEEK_SET) == 0;
#endif
  }

private:
  fileptr_t f;
};

struct MemoryStream final : public AnyStream {
  MemoryStream(const char* start_, size_t size)
    : start(start_), end(start_ + size), cur(start_) {}

  char* gets(char* line, int size) override {
    --size; // fgets reads in at most one less than size characters
    if (cur >= end)
      return nullptr;
    if (size > end - cur)
      size = int(end - cur);
    const char* nl = (const char*) std::memchr(cur, '\n', size);
    size_t len = nl ? nl - cur + 1 : size;
    std::memcpy(line, cur, len);
    line[len] = '\0';
    cur += len;
    return line;
  }
  int getc() override { return cur < end ? *cur++ : EOF; }

  bool read(void* buf, size_t len) override {
    if (cur + len > end)
      return false;
    std::memcpy(buf, cur, len);
    cur += len;
    return true;
  }

  std::string read_rest() override {
    const char* last = cur;
    cur = end;
    return std::string(last, end);
  }

  bool seek(std::ptrdiff_t offset) override {
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

  const std::string& path() const { return path_; }
  const std::string& basepath() const { return path_; }

  // Does the path stands for stdin?
  // Each reading function needs to call it (some functions use stdin
  // and some std::cin, so we don't try to unify it here).
  bool is_stdin() const { return path() == "-"; }

  // providing the same interface as MaybeGzipped
  bool is_compressed() const { return false; }
  // for reading (uncompressing into memory) the whole file at once
  CharArray uncompress_into_buffer(size_t=0) { return {}; }

  std::unique_ptr<AnyStream> create_stream() {
    return std::unique_ptr<AnyStream>(new FileStream(path().c_str(), "rb"));
  }

private:
  std::string path_;
};

} // namespace gemmi
#endif
