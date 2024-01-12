// Copyright 2017 Global Phasing Ltd.
//
// Functions for transparent reading of gzipped files. Uses zlib.

#ifndef GEMMI_GZ_HPP_
#define GEMMI_GZ_HPP_
#include <string>
#include "fail.hpp"     // GEMMI_DLL
#include "input.hpp"    // BasicInput
#include "util.hpp"     // iends_with

namespace gemmi {

GEMMI_DLL size_t estimate_uncompressed_size(const std::string& path);

// the same interface as FileStream and MemoryStream
struct GEMMI_DLL GzStream {
  void* f;  // implementation detail
  char* gets(char* line, int size);
  int getc();
  bool read(void* buf, size_t len);
};

class GEMMI_DLL MaybeGzipped : public BasicInput {
public:
  explicit MaybeGzipped(const std::string& path);
  ~MaybeGzipped();
  size_t gzread_checked(void* buf, size_t len);
  bool is_compressed() const { return iends_with(path(), ".gz"); }
  std::string basepath() const {
    return is_compressed() ? path().substr(0, path().size() - 3) : path();
  }

  CharArray uncompress_into_buffer(size_t limit=0);
  GzStream get_uncompressing_stream();

private:
  void* file_ = nullptr;
};

} // namespace gemmi

#endif
