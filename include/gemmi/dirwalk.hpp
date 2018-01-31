// Copyright 2018 Global Phasing Ltd.
//
// Two class for iterating files in a directory tree, top-down,
// in an alphabetical order.  It wraps the tinydir library (as we cannot
// depend on C++17 <filesystem> yet).
// DirWalk iterates through all files and directories.
// CifWalk yields only cif files (either files that end with .cif or .cif.gz,
// or files that look like SF mmCIF files from wwPDB, e.g. r3aaasf.ent.gz).
// It's good for traversing a local copy of the wwPDB archive.
//
// Usage:
//   for (tinydir_file& file : gemmi::DirWalk(top_dir))
//     do_something(file);
// or
//   for (const char* file : gemmi::CifWalk(top_dir))
//     do_something(file);
// You should also catch std::runtime_error.

#ifndef GEMMI_DIRWALK_HPP_
#define GEMMI_DIRWALK_HPP_

#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>
#include <cassert>

#if defined(_MSC_VER) && !defined(NOMINMAX)
# define NOMINMAX
#endif
#include <tinydir.h>

namespace gemmi {

class DirWalk {
public:
  explicit DirWalk(const char* path) {
    if (tinydir_file_open(&top_, path) == -1) {
      //std::perror(nullptr);
      throw std::runtime_error("Cannot open file or directory: " +
                               std::string(path));
    }
  }
  explicit DirWalk(const std::string& path) { DirWalk(path.c_str()); }
  ~DirWalk() {
    for (auto& d : dirs_)
      tinydir_close(&d.second);
  }
  void push_dir(int cur_pos, const char* path) {
    dirs_.emplace_back();
    dirs_.back().first = cur_pos;
    if (tinydir_open_sorted(&dirs_.back().second, path) == -1)
      throw std::runtime_error("Cannot open directory: " + std::string(path));
  }
  int pop_dir() {
    assert(!dirs_.empty());
    int old_pos = dirs_.back().first;
    tinydir_close(&dirs_.back().second);
    dirs_.pop_back();
    return old_pos;
  }

  struct Iter {
    DirWalk& walk;
    size_t cur;

    const tinydir_dir& get_dir() const { return walk.dirs_.back().second; }

    const tinydir_file& operator*() const {
      if (walk.dirs_.empty())
        return walk.top_;
      assert(cur < get_dir().n_files);
      return get_dir()._files[cur];
    }

    bool is_special(const char* name) const {
      return strcmp(name, ".") == 0 || strcmp(name, "..") == 0;
    }

    size_t depth() const { return walk.dirs_.size(); }

    void operator++() { // depth first
      const tinydir_file& tf = **this;
      if (tf.is_dir) {
        walk.push_dir(cur, tf.path);
        cur = 0;
      } else {
        cur++;
      }
      while (!walk.dirs_.empty()) {
        if (cur == get_dir().n_files)
          cur = walk.pop_dir() + 1;
        else if (is_special(get_dir()._files[cur].name))
          cur++;
        else
          break;
      }
    }

    bool operator!=(const Iter& o) const {  // it's only for "!= end()"
      return !(depth() == 0 && cur == o.cur);
    }
  };

  Iter begin() { return Iter{*this, 0}; }
  Iter end() { return Iter{*this, 1}; }
  bool is_single_file() { return !top_.is_dir; }

private:
  friend struct Iter;
  tinydir_file top_;
  std::vector<std::pair<int, tinydir_dir>> dirs_;
};

struct IsCifFile {
  static bool is_ok(const std::string& filename) {
    return gemmi::giends_with(filename, ".cif") ||
        // the SF mmCIF files from PDB don't have the "cif" extension,
        // they have names such as divided/structure_factors/aa/r3aaasf.ent.gz
        (filename[0] == 'r' && gemmi::giends_with(filename, "sf.ent"));
  }
};

struct IsPdbFile {
  static bool is_ok(const std::string& filename) {
    return gemmi::giends_with(filename, ".pdb") ||
           (gemmi::giends_with(filename, ".ent") &&
            !(filename[0] == 'r' && gemmi::giends_with(filename, "sf.ent")));
  }
};


template<typename Check>
class FileWalk : public DirWalk {
public:
  explicit FileWalk(const char* path) : DirWalk(path) {}
  explicit FileWalk(const std::string& path) : FileWalk(path.c_str()) {}

  struct CifIter : Iter {
    CifIter(Iter&& iter) : Iter(iter) {}
    void operator++() {
      for (;;) {
        Iter::operator++();
        const tinydir_file& f = Iter::operator*();
        if ((!f.is_dir && Check::is_ok(f.name))
            || walk.is_single_file()
            || (depth() == 0 && cur == 1))
          break;
      }
    }
    const char* operator*() const { return Iter::operator*().path; }
  };
  CifIter begin() {
    CifIter it = DirWalk::begin();
    if (!is_single_file()) // i.e. the top item is a directory
      ++it;
    return it;
  }
  CifIter end() { return DirWalk::end(); }
};

using CifWalk = FileWalk<IsCifFile>;
using PdbWalk = FileWalk<IsPdbFile>;

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
