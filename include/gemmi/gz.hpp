//! @file
//! @brief Functions for transparent reading of gzipped files.
//!
//! Functions for transparent reading of gzipped files. Uses zlib.

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

//! Description of linked zlib library.
GEMMI_DLL extern const char* const zlib_description;

//! @brief Estimate uncompressed size of gzipped file.
//! @param path Path to gzipped file
//! @return Estimated uncompressed size (may be inaccurate)
GEMMI_DLL size_t estimate_uncompressed_size(const std::string& path);

//! @brief Stream for reading gzipped data.
//!
//! The same interface as FileStream and MemoryStream.
struct GEMMI_DLL GzStream final : public AnyStream {
  GzStream(void* f_) : f(f_) {}
  char* gets(char* line, int size) override;
  int getc() override;
  bool read(void* buf, size_t len) override;
  bool skip(size_t n) override;
  long tell() override;
  std::string read_rest() override;

private:
  void* f;  //!< Implementation detail (gzFile)
};

//! @brief Input that may be gzipped or plain.
//!
//! Handles both compressed (.gz) and uncompressed files transparently.
class GEMMI_DLL MaybeGzipped : public BasicInput {
public:
  //! @brief Construct from file path.
  //! @param path Path to file (gzipped or plain)
  explicit MaybeGzipped(const std::string& path);

  ~MaybeGzipped();

  //! @brief Read and decompress data.
  //! @param buf Output buffer
  //! @param len Number of bytes to read
  //! @return Number of bytes read
  size_t gzread_checked(void* buf, size_t len);

  //! @brief Check if file is gzipped.
  //! @return True if path ends with .gz
  bool is_compressed() const { return iends_with(path(), ".gz"); }

  //! @brief Get path without .gz extension.
  //! @return Base path (without .gz if present)
  std::string basepath() const {
    return is_compressed() ? path().substr(0, path().size() - 3) : path();
  }

  //! @brief Decompress entire file into memory buffer.
  //! @param limit Optional size limit (0 = no limit)
  //! @return Buffer containing uncompressed data
  CharArray uncompress_into_buffer(size_t limit=0);

  //! @brief Create stream for reading file.
  //! @return Stream (GzStream for .gz, FileStream otherwise)
  std::unique_ptr<AnyStream> create_stream();

private:
  void* file_ = nullptr;  //!< Gzip file handle (if compressed)
};

} // namespace gemmi

#endif
