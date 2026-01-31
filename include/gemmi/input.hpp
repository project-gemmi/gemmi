//! @file
//! @brief Input abstraction for file reading and decompression.
//!
//! Input abstraction.
//! Used to decouple file reading and decompression.

// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.
// Used to decouple file reading and decompression.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <cstddef> // for ptrdiff_t
#include <cstdio>  // for FILE, fseek, fread
#include <cstring> // for memchr
#include <string>
#include "fileutil.hpp"  // for fileptr_t

namespace gemmi {

//! @brief Base class for FileStream, MemoryStream and GzStream.
//!
//! Provides uniform interface for reading from files, memory, or gzipped streams.
struct AnyStream {
  virtual ~AnyStream() = default;

  //! @brief Read line into buffer (for pdb, copy_line()).
  //! @param line Output buffer
  //! @param size Buffer size
  //! @return Pointer to line or nullptr on EOF
  virtual char* gets(char* line, int size) = 0;

  //! @brief Read single character (for copy_line()).
  //! @return Character or EOF
  virtual int getc() = 0;

  //! @brief Read binary data (for ccp4, mtz).
  //! @param buf Output buffer
  //! @param len Number of bytes to read
  //! @return True if successfully read len bytes
  virtual bool read(void* buf, size_t len) = 0;

  //! @brief Get current position (temporary, for testing).
  //! @return Current position
  //!
  //! These are not used in GzStream because MemoryStream is used for mtz.
  virtual long tell() = 0;

  //! @brief Skip bytes (for reading mtz without data).
  //! @param n Number of bytes to skip
  //! @return True on success
  virtual bool skip(size_t n) = 0;

  //! @brief Read remaining data (for mtz appendix).
  //! @return Remaining data as string
  virtual std::string read_rest() { return {}; }

  //! @brief Copy line handling overflow (for pdb, xds_ascii).
  //! @param line Output buffer
  //! @param size Buffer size
  //! @return Number of characters read
  //!
  //! If a line is longer than size we discard the rest of it.
  size_t copy_line(char* line, int size) {
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

//! @brief Stream reading from FILE*.
struct FileStream final : public AnyStream {
  //! @brief Construct from existing FILE*.
  //! @param f_ File pointer (not owned)
  FileStream(std::FILE* f_) : f(f_, needs_fclose{false}) {}

  //! @brief Construct from path.
  //! @param path File path (or "-" for stdin if allowed)
  //! @param mode File mode ("rb", "r", etc.)
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

  long tell() override {
    return std::ftell(f.get());
  }

  bool skip(size_t n) override {
#if defined(_MSC_VER)
    int result = _fseeki64(f.get(), (std::ptrdiff_t)n, SEEK_CUR);
#elif defined(__MINGW32__)
    int result = fseeko(f.get(), (_off_t)n, SEEK_CUR);
#else
    int result = std::fseek(f.get(), (long)n, SEEK_CUR);
#endif
    if (result != 0) {
      char buf[512];
      while (n >= sizeof(buf)) {
        if (std::fread(buf, sizeof(buf), 1, f.get()) != 1)
          return false;
        n -= sizeof(buf);
      }
      if (n > 0 && std::fread(buf, n, 1, f.get()) != 1)
        return false;
    }
    return true;
  }

private:
  fileptr_t f;  //!< File pointer
};

//! @brief Stream reading from memory buffer.
struct MemoryStream final : public AnyStream {
  //! @brief Construct from memory buffer.
  //! @param start_ Pointer to start of buffer
  //! @param size Buffer size
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

  long tell() override {
    return cur - start;
  }
  bool skip(size_t n) override {
    cur += n;
    return cur < end;
  }

private:
  const char* const start;  //!< Start of buffer
  const char* const end;    //!< End of buffer
  const char* cur;          //!< Current position
};

//! @brief Basic input from file path.
//!
//! Provides the same interface as MaybeGzipped for non-compressed files.
class BasicInput {
public:
  //! @brief Construct from file path.
  //! @param path File path (or "-" for stdin)
  explicit BasicInput(const std::string& path) : path_(path) {}

  //! @brief Get file path.
  //! @return File path
  const std::string& path() const { return path_; }

  //! @brief Get base path (same as path() for non-compressed).
  //! @return Base path
  const std::string& basepath() const { return path_; }

  //! @brief Check if path represents stdin.
  //! @return True if path is "-"
  //!
  //! Does the path stands for stdin?
  //! Each reading function needs to call it (some functions use stdin
  //! and some std::cin, so we don't try to unify it here).
  bool is_stdin() const { return path() == "-"; }

  //! @brief Check if file is compressed.
  //! @return False (providing the same interface as MaybeGzipped)
  bool is_compressed() const { return false; }

  //! @brief Uncompress file into buffer (no-op for BasicInput).
  //! @return Empty buffer (for reading uncompressing into memory the whole file at once)
  CharArray uncompress_into_buffer(size_t=0) { return {}; }

  //! @brief Create stream for reading file.
  //! @return FileStream for the path
  std::unique_ptr<AnyStream> create_stream() {
    return std::unique_ptr<AnyStream>(new FileStream(path().c_str(), "rb"));
  }

private:
  std::string path_;  //!< File path
};

} // namespace gemmi
#endif
