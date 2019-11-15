// Copyright 2019 Global Phasing Ltd.
//
// Ofstream - wrapper around std::ofstream with two extra features:
//  - on MSVC supports Unicode filenames (the filename is passed in UTF-8)
//  - Ofstream("-", ...) returns wrapper around the stream in the 2nd arg
//    (which is used to interpret "-" as stdout or stderr).

#ifndef GEMMI_OFSTREAM_HPP_
#define GEMMI_OFSTREAM_HPP_

#if defined(_MSC_VER) && !defined(GEMMI_USE_FOPEN)
# include "utf.hpp"
#endif
#include <fstream>
#include <memory>
#include "fail.hpp"

namespace gemmi {

// note: move of std::ofstream doesn't work in GCC 4.8.

struct Ofstream {
  Ofstream(const std::string& filename, std::ostream* dash=nullptr) {
    if (filename.size() == 1 && filename[0] == '-' && dash) {
      ptr_ = dash;
      return;
    }
    keeper_.reset(new std::ofstream);
#if defined(_MSC_VER) && !defined(GEMMI_USE_FOPEN)
    std::wstring wfilename = UTF8_to_wchar(filename.c_str());
    keeper_->open(wfilename.c_str());
#else
    keeper_->open(filename);
#endif
    if (!*keeper_)
      fail("Failed to open file for writing: " + filename);
    ptr_ = keeper_.get();
  }

  std::ostream* operator->() { return ptr_; }
  std::ostream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ofstream> keeper_;
  std::ostream* ptr_;
};


} // namespace gemmi
#endif
