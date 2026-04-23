// Copyright 2019 Global Phasing Ltd.
//
// Ofstream and Ifstream: wrappers around std::ofstream and std::ifstream.
//
// They offer two extra features:
//  - on MSVC supports Unicode filenames (the filename is passed in UTF-8),
//  - optionally, filename "-" can be interpreted as stdout or stderr.

#ifndef GEMMI_OFSTREAM_HPP_
#define GEMMI_OFSTREAM_HPP_

#if defined(_MSC_VER)
# include "utf.hpp"
#elif defined(_WIN32) && defined(__has_include)
// including <filesystem> in MinGW 8 gives error (in bits/fs_path.h:237:47),
// which was fixed in GCC 9
# if __has_include(<filesystem>) && !(defined(__MINGW32__) && __GNUC__ < 9)
#  include <filesystem>
#  include "utf.hpp"
# endif
#endif
#include <fstream>
#include <memory>
#include "fail.hpp"

namespace gemmi {

/// Open a file stream with UTF-8 filename support.
/// @brief Helper to open streams with UTF-8 paths on Windows.
/// @tparam T stream type (std::ofstream or std::ifstream)
/// @param ptr pointer to stream object to open
/// @param filename UTF-8 encoded filename
template<typename T>
inline void open_stream_from_utf8_path(T& ptr, const std::string& filename) {
#if defined(_MSC_VER)
    std::wstring wfilename = UTF8_to_wchar(filename.c_str());
    ptr->open(wfilename.c_str());
#elif defined(_WIN32) && defined(__cpp_lib_filesystem)
    std::wstring wfilename = UTF8_to_wchar(filename.c_str());
    ptr->open(std::filesystem::path(wfilename));
#else
    ptr->open(filename);
#endif
}

// note: move of std::ofstream doesn't work in GCC 4.8.

/// @brief Output file stream wrapper with UTF-8 filename support.
/// @details Handles filename "-" as stdout and UTF-8 paths on Windows.
struct Ofstream {
  /// Open output file with optional dash handling.
  /// @brief Construct output stream.
  /// @param filename UTF-8 file path (or "-" to use dash_stream)
  /// @param dash pointer to stream to use if filename is "-" (typically std::cout)
  Ofstream(const std::string& filename, std::ostream* dash=nullptr) {
    if (filename.size() == 1 && filename[0] == '-' && dash) {
      ptr_ = dash;
      return;
    }
    keeper_.reset(new std::ofstream);
    open_stream_from_utf8_path(keeper_, filename);
    if (!*keeper_)
      sys_fail("Failed to open " + filename + " for writing");
    ptr_ = keeper_.get();
  }

  /// Get pointer to stream.
  /// @brief Access stream object via pointer.
  /// @return pointer to std::ostream
  std::ostream* operator->() { return ptr_; }
  /// Get reference to stream.
  /// @brief Access stream object by reference.
  /// @return reference to std::ostream
  std::ostream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ofstream> keeper_;
  std::ostream* ptr_;
};

/// @brief Input file stream wrapper with UTF-8 filename support.
/// @details Handles filename "-" as stdin and UTF-8 paths on Windows.
struct Ifstream {
  /// Open input file with optional dash handling.
  /// @brief Construct input stream.
  /// @param filename UTF-8 file path (or "-" to use dash_stream)
  /// @param dash pointer to stream to use if filename is "-" (typically std::cin)
  Ifstream(const std::string& filename, std::istream* dash=nullptr) {
    if (filename.size() == 1 && filename[0] == '-' && dash) {
      ptr_ = dash;
      return;
    }
    keeper_.reset(new std::ifstream);
    open_stream_from_utf8_path(keeper_, filename);
    if (!*keeper_)
      sys_fail("Failed to open " + filename);
    ptr_ = keeper_.get();
  }

  /// Get pointer to stream.
  /// @brief Access stream object via pointer.
  /// @return pointer to std::istream
  std::istream* operator->() { return ptr_; }
  /// Get reference to stream.
  /// @brief Access stream object by reference.
  /// @return reference to std::istream
  std::istream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ifstream> keeper_;
  std::istream* ptr_;
};



} // namespace gemmi
#endif
