// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <memory>  // for unique_ptr
#include <string>

namespace gemmi {

class JustFile {
public:
  struct DummyLineInput {
    explicit operator bool() const { return false; }
    char* gets(char*, int) { return nullptr; }
    int getc() { return 0; }
  };
  explicit JustFile(const std::string& path) : path_(path) {}
  bool is_stdin() const { return false; };
  const std::string& path() const { return path_; };
  std::unique_ptr<char[]> memory() { return nullptr; }
  size_t memory_size() const { return 0; };
  DummyLineInput get_line_stream() const { return DummyLineInput(); }
  //bool get_line_stream() const { return false; }
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
// vim:sw=2:ts=2:et
