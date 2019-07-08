// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.
// Used to decouple file reading and uncompression.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <memory>  // for unique_ptr
#include <string>

namespace gemmi {

class JustFile {
public:
  // for incremental reading, used with get_stream()
  struct DummyStream {
    explicit operator bool() const { return false; }
    char* gets(char*, int) { return nullptr; }  // used in pdb.hpp
    int getc() { return 0; }                    // used in pdb.hpp
    bool read(void*, int) { return false; }     // used in ccp4.hpp
  };

  explicit JustFile(const std::string& path) : path_(path) {}

  const std::string& path() const { return path_; };

  // for handling of some file names as stdin, ("-" in MaybeStdin);
  // each reading function needs to call it (some functions use stdin
  // and some std::cin, so we don't try to unify it here)
  bool is_stdin() const { return false; };

  // for reading (uncompressing into memory) the whole file at once
  std::unique_ptr<char[]> memory() { return nullptr; }
  size_t memory_size() const { return 0; };

  // for incremental reading
  DummyStream get_stream() const { return DummyStream(); }

private:
  std::string path_;
};

// for transparent handling of stdin along filenames
class MaybeStdin : public JustFile {
public:
  explicit MaybeStdin(const std::string& path) : JustFile(path) {}
  bool is_stdin() const { return path() == "-"; };
};

} // namespace gemmi
#endif
