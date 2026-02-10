//! @file
//! @brief Stream wrappers with Unicode and dash support.
//!
//! Ofstream and Ifstream: wrappers around std::ofstream and std::ifstream.
//!
//! They offer two extra features:
//!  - on MSVC supports Unicode filenames (the filename is passed in UTF-8),
//!  - optionally, filename "-" can be interpreted as stdout or stderr.

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

//! @brief Open stream from UTF-8 path.
//! @tparam T Smart pointer type to stream
//! @param ptr Smart pointer to stream
//! @param filename UTF-8 encoded filename
//!
//! Helper function that handles Unicode paths on Windows.
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

//! @brief Output file stream wrapper with Unicode and dash support.
//!
//! Wraps std::ofstream with UTF-8 path support on Windows and optional
//! interpretation of "-" as stdout/stderr.
struct Ofstream {
  //! @brief Construct output stream.
  //! @param filename Filename (or "-" for dash stream)
  //! @param dash Optional stream to use for "-" (e.g., stdout or stderr)
  //! @throws std::system_error if file cannot be opened
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

  //! @brief Get pointer to underlying stream.
  //! @return Stream pointer
  std::ostream* operator->() { return ptr_; }

  //! @brief Get reference to underlying stream.
  //! @return Stream reference
  std::ostream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ofstream> keeper_;  //!< Owned stream
  std::ostream* ptr_;                      //!< Pointer to stream
};

//! @brief Input file stream wrapper with Unicode and dash support.
//!
//! Wraps std::ifstream with UTF-8 path support on Windows and optional
//! interpretation of "-" as stdin.
struct Ifstream {
  //! @brief Construct input stream.
  //! @param filename Filename (or "-" for dash stream)
  //! @param dash Optional stream to use for "-" (e.g., stdin)
  //! @throws std::system_error if file cannot be opened
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

  //! @brief Get pointer to underlying stream.
  //! @return Stream pointer
  std::istream* operator->() { return ptr_; }

  //! @brief Get reference to underlying stream.
  //! @return Stream reference
  std::istream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ifstream> keeper_;  //!< Owned stream
  std::istream* ptr_;                      //!< Pointer to stream
};



} // namespace gemmi
#endif
