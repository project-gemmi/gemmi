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

GEMMI_DLL extern const char* const zlib_description;

GEMMI_DLL size_t estimate_uncompressed_size(const std::string& path);

// the same interface as FileStream and MemoryStream
struct GEMMI_DLL GzStream final : public AnyStream {
  GzStream(void* f_) : f(f_) {}
  char* gets(char* line, int size) override;
  int getc() override;
  bool read(void* buf, size_t len) override;
private:
  void* f;  // implementation detail
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

  std::unique_ptr<AnyStream> create_stream();

private:
  void* file_ = nullptr;
};

} // namespace gemmi

#endif
