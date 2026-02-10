//! @file
//! @brief File-related utilities for reading and manipulating files.
//!
//! Provides functions for file operations including opening files (with Windows Unicode
//! support), determining file sizes, byte swapping for endianness, and reading files
//! into memory buffers. Includes CharArray class for managing dynamically allocated buffers.

// Copyright 2018 Global Phasing Ltd.

#ifndef GEMMI_FILEUTIL_HPP_
#define GEMMI_FILEUTIL_HPP_

#include <cassert>
#include <cstdio>    // for FILE, fopen, fclose
#include <cstdint>
#include <cstdlib>   // for malloc, realloc
#include <cstring>   // for strlen
#include <initializer_list>
#include <memory>    // for unique_ptr
#include "fail.hpp"  // for sys_fail

#if defined(_WIN32) && !defined(GEMMI_USE_FOPEN)
#include "utf.hpp"
#endif

namespace gemmi {

//! @brief Extract basename from path and strip directory and specified suffixes.
//! @param path File path
//! @param exts List of extensions to remove
//! @return Basename without directory or specified extensions
//!
//! Strip directory and suffixes from filename.
inline std::string path_basename(const std::string& path,
                                 std::initializer_list<const char*> exts) {
  size_t pos = path.find_last_of("\\/");
  std::string basename = pos == std::string::npos ? path : path.substr(pos + 1);
  for (const char* ext : exts) {
    size_t len = std::strlen(ext);
    if (basename.size() > len &&
        basename.compare(basename.length() - len, len, ext, len) == 0)
      basename.resize(basename.length() - len);
  }
  return basename;
}

// file operations

//! @brief Custom deleter for FILE* in unique_ptr.
//!
//! Deleter for fileptr_t.
struct needs_fclose {
  bool use_fclose;  //!< Whether to call fclose (false for stdin/stdout)
  void operator()(std::FILE* f) const noexcept {
    if (use_fclose)
      std::fclose(f);
  }
};

//! Unique pointer to FILE with custom deleter.
typedef std::unique_ptr<std::FILE, needs_fclose> fileptr_t;

//! @brief Open file with error handling.
//! @param path File path
//! @param mode Open mode (as for fopen)
//! @return Unique pointer to opened FILE
//! @throws std::system_error if file cannot be opened
//!
//! On Windows, uses wide character path for Unicode support.
inline fileptr_t file_open(const char* path, const char* mode) {
  std::FILE* file;
#if defined(_WIN32) && !defined(GEMMI_USE_FOPEN)
  std::wstring wpath = UTF8_to_wchar(path);
  std::wstring wmode = UTF8_to_wchar(mode);
  if ((file = ::_wfopen(wpath.c_str(), wmode.c_str())) == nullptr)
#else
  if ((file = std::fopen(path, mode)) == nullptr)
#endif
    sys_fail(std::string("Failed to open ") + path +
             (*mode == 'w' ? " for writing" : ""));
  return fileptr_t(file, needs_fclose{true});
}

//! @brief Open file or use stdin/stdout for "-".
//! @param path File path (or "-" for stdin/stdout)
//! @param mode Open mode
//! @param dash_stream Stream to use for "-" (stdin or stdout)
//! @return Unique pointer to opened FILE
//!
//! Helper function for treating "-" as stdin or stdout.
inline fileptr_t file_open_or(const char* path, const char* mode,
                              std::FILE* dash_stream) {
  if (path[0] == '-' && path[1] == '\0')
    return fileptr_t(dash_stream, needs_fclose{false});
  return file_open(path, mode);
}

//! @brief Determine file size.
//! @param f Open file handle
//! @param path File path (for error messages)
//! @return File size in bytes
//! @throws std::system_error if seek/tell operations fail
inline std::size_t file_size(std::FILE* f, const std::string& path) {
  if (std::fseek(f, 0, SEEK_END) != 0)
    sys_fail(path + ": fseek failed");
  long length = std::ftell(f);
  if (length < 0)
    sys_fail(path + ": ftell failed");
  if (std::fseek(f, 0, SEEK_SET) != 0)
    sys_fail(path + ": fseek failed");
  return length;
}

//! @brief Check if system is little-endian.
//! @return True if little-endian
//!
//! Helper function for working with binary files.
inline bool is_little_endian() {
  std::uint32_t x = 1;
  return *reinterpret_cast<char *>(&x) == 1;
}

//! @brief Swap byte order of 2-byte value.
//! @param start Pointer to value
inline void swap_two_bytes(void* start) {
  char* bytes = static_cast<char*>(start);
  std::swap(bytes[0], bytes[1]);
}

//! @brief Swap byte order of 4-byte value.
//! @param start Pointer to value
inline void swap_four_bytes(void* start) {
  char* bytes = static_cast<char*>(start);
  std::swap(bytes[0], bytes[3]);
  std::swap(bytes[1], bytes[2]);
}

//! @brief Swap byte order of 8-byte value.
//! @param start Pointer to value
inline void swap_eight_bytes(void* start) {
  char* bytes = static_cast<char*>(start);
  std::swap(bytes[0], bytes[7]);
  std::swap(bytes[1], bytes[6]);
  std::swap(bytes[2], bytes[5]);
  std::swap(bytes[3], bytes[4]);
}

//! @brief Dynamically allocated character array with automatic memory management.
//!
//! Wraps malloc/realloc/free for C-style character buffers. Uses unique_ptr
//! with custom deleter for automatic cleanup.
class CharArray {
  std::unique_ptr<char, decltype(&std::free)> ptr_;
  size_t size_;
public:
  //! @brief Construct empty array.
  CharArray() : ptr_(nullptr, &std::free), size_(0) {}

  //! @brief Construct array of specified size.
  //! @param n Size in bytes
  explicit CharArray(size_t n) : ptr_((char*)std::malloc(n), &std::free), size_(n) {}

  //! @brief Check if array is allocated.
  //! @return True if non-null
  explicit operator bool() const { return (bool)ptr_; }

  //! @brief Get pointer to data.
  //! @return Mutable pointer
  char* data() { return ptr_.get(); }

  //! @brief Get const pointer to data.
  //! @return Const pointer
  const char* data() const { return ptr_.get(); }

  //! @brief Get size of array.
  //! @return Size in bytes
  size_t size() const { return size_; }

  //! @brief Set size (without reallocation).
  //! @param n New size
  void set_size(size_t n) { size_ = n; }

  //! @brief Resize array (may reallocate).
  //! @param n New size in bytes
  //! @throws std::runtime_error if out of memory
  void resize(size_t n) {
    char* new_ptr = (char*) std::realloc(ptr_.get(), n);
    if (!new_ptr && n != 0)
      fail("Out of memory.");
    (void) ptr_.release();  // NOLINT(bugprone-unused-return-value)
    ptr_.reset(new_ptr);
    size_ = n;
  }

  //! @brief Remove first n bytes making space for more text at the returned position.
  //! @param n Number of bytes to remove
  //! @return Pointer to position after kept data
  char* roll(size_t n) {
    assert(n <= size());
    std::memmove(data(), data() + n, n);
    return data() + n;
  }
};


//! @brief Read file into memory buffer.
//! @param path File path
//! @return CharArray containing file contents
//! @throws std::system_error if file operations fail
//!
//! Reading file into a memory buffer (optimized: uses fseek to determine file size).
inline CharArray read_file_into_buffer(const std::string& path) {
  fileptr_t f = file_open(path.c_str(), "rb");
  size_t size = file_size(f.get(), path);
  CharArray buffer(size);
  if (std::fread(buffer.data(), size, 1, f.get()) != 1)
    sys_fail(path + ": fread failed");
  return buffer;
}

//! @brief Read stdin into memory buffer.
//! @return CharArray containing stdin contents
inline CharArray read_stdin_into_buffer() {
  size_t n = 0;
  CharArray buffer(16 * 1024);
  for (;;) {
    n += std::fread(buffer.data() + n, 1, buffer.size() - n, stdin);
    if (n != buffer.size()) {
      buffer.set_size(n);
      break;
    }
    buffer.resize(2*n);
  }
  return buffer;
}

//! @brief Read input into buffer (file, stdin, or compressed).
//! @tparam T Input type with is_compressed(), is_stdin(), path() methods
//! @param input Input object
//! @return CharArray containing input contents
template<typename T>
inline CharArray read_into_buffer(T&& input) {
  if (input.is_compressed())
    return input.uncompress_into_buffer();
  if (input.is_stdin())
    return read_stdin_into_buffer();
  return read_file_into_buffer(input.path());
}

} // namespace gemmi
#endif
