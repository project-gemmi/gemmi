//! @file
//! @brief Classes for iterating over files in a directory tree.
//!
//! Classes for iterating over files in a directory tree, top-down,
//! in alphabetical order. Wraps the tinydir library (as we cannot yet
//! depend on C++17 <filesystem>).
//!
//! DirWalk<> iterates through all files and directories.
//! CifWalk yields only cif files (either files that end with .cif or .cif.gz,
//! or files that look like SF mmCIF files from wwPDB, e.g. r3aaasf.ent.gz).
//! It's good for traversing a local copy of the wwPDB archive.
//! PdbWalk: .pdb or .ent (optionally with .gz) except r????sf.ent
//! CoorFileWalk: .cif, .pdb or .ent (optionally with .gz)
//!               except r????sf.ent and *-sf.cif
//!
//! Usage:
//!   for (const std::string& file : gemmi::DirWalk<>(top_dir))
//!     do_something(file);
//! or
//!   for (const std::string& file : gemmi::CifWalk(top_dir))
//!     do_something(file);
//! You should also catch std::runtime_error.

// Copyright 2018 Global Phasing Ltd.
//
// Classes for iterating over files in a directory tree, top-down,
// in alphabetical order. Wraps the tinydir library (as we cannot yet
// depend on C++17 <filesystem>).

// DirWalk<> iterates through all files and directories.
// CifWalk yields only cif files (either files that end with .cif or .cif.gz,
// or files that look like SF mmCIF files from wwPDB, e.g. r3aaasf.ent.gz).
// It's good for traversing a local copy of the wwPDB archive.
// PdbWalk: .pdb or .ent (optionally with .gz) except r????sf.ent
// CoorFileWalk: .cif, .pdb or .ent (optionally with .gz)
//               except r????sf.ent and *-sf.cif
//
// Usage:
//   for (const std::string& file : gemmi::DirWalk<>(top_dir))
//     do_something(file);
// or
//   for (const std::string& file : gemmi::CifWalk(top_dir))
//     do_something(file);
// You should also catch std::runtime_error.

#ifndef GEMMI_DIRWALK_HPP_
#define GEMMI_DIRWALK_HPP_

#include <string>
#include <vector>
#include <cassert>
#if defined(_MSC_VER) && !defined(NOMINMAX)
# define NOMINMAX
#endif
#include "third_party/tinydir.h"

#include "util.hpp"  // for giends_with
#include "fail.hpp"  // for sys_fail
#include "pdb_id.hpp" // for is_pdb_code, expand_pdb_code_to_path
#include "glob.hpp"  // for glob_match
#if defined(_WIN32) && defined(_UNICODE)
 #include "utf.hpp"
#endif

namespace gemmi {

//! @brief Convert platform path to UTF-8 string.
//! @param path Platform-specific path
//! @return UTF-8 encoded path
inline std::string as_utf8(const _tinydir_char_t* path) {
#if defined(_WIN32) && defined(_UNICODE)
  return wchar_to_UTF8(path);
#else
  return path;
#endif
}


namespace impl {

//! @brief Check if filename matches SF mmCIF pattern from PDB.
//! @param filename File name to check
//! @return True if matches pattern like r3aaasf.ent
//!
//! The SF mmCIF files from PDB have names such as
//! divided/structure_factors/aa/r3aaasf.ent.gz
inline bool is_rxsf_ent_filename(const std::string& filename) {
  return filename[0] == 'r' && giends_with(filename, "sf.ent")
         && filename.find('.') >= 4;
}

//! @brief Filter for mmCIF files.
//!
//! Actually we don't know what kind of cif file it is.
struct IsMmCifFile {
  static bool check(const std::string& filename) {
    return giends_with(filename, ".cif") || giends_with(filename, ".mmcif");
  }
};

//! @brief Filter for CIF files (including SF mmCIF).
struct IsCifFile {
  static bool check(const std::string& filename) {
    return giends_with(filename, ".cif") || is_rxsf_ent_filename(filename);
  }
};

//! @brief Filter for PDB files.
struct IsPdbFile {
  static bool check(const std::string& filename) {
    return giends_with(filename, ".pdb") ||
           (giends_with(filename, ".ent") && !is_rxsf_ent_filename(filename));
  }
};

//! @brief Filter for coordinate files (excluding SF files).
struct IsCoordinateFile {
  static bool check(const std::string& filename) {
    // the SF mmCIF files from RCSB website have names such as 3AAA-sf.cif
    return IsPdbFile::check(filename) ||
           (IsMmCifFile::check(filename) && !giends_with(filename, "-sf.cif"));
  }
};

//! @brief Filter that accepts any file.
struct IsAnyFile {
  static bool check(const std::string&) { return true; }
};

//! @brief Filter using glob pattern.
struct IsMatchingFile {
  bool check(const std::string& filename) const {
    return glob_match(pattern, filename);
  }
  std::string pattern;  //!< Glob pattern
};

//! @brief Open tinydir file with UTF-8 path.
//! @param file Tinydir file structure
//! @param path UTF-8 encoded path
//! @return 0 on success, -1 on failure
inline int utf8_tinydir_file_open(tinydir_file* file, const char* path) {
#if defined(_WIN32) && defined(_UNICODE)
  return tinydir_file_open(file, UTF8_to_wchar(path).c_str());
#else
  return tinydir_file_open(file, path);
#endif
}

} // namespace impl


//! @brief Directory tree iterator with filtering.
//! @tparam FileOnly If true, iterate only files (not directories)
//! @tparam Filter File filter predicate
//!
//! Iterates depth-first through directory tree in alphabetical order.
template<bool FileOnly=true, typename Filter=impl::IsAnyFile>
class DirWalk {
public:
  //! @brief Construct directory walker.
  //! @param path Directory path or PDB code
  //! @param try_pdbid Character for PDB code expansion ('M', 'S', etc.)
  //! @throws std::system_error if path cannot be opened
  explicit DirWalk(const char* path, char try_pdbid='\0') {
    if (impl::utf8_tinydir_file_open(&top_, path) != -1)
      return;
    if (try_pdbid != '\0' && is_pdb_code(path)) {
      std::string epath = expand_pdb_code_to_path(path, try_pdbid, true);
      if (impl::utf8_tinydir_file_open(&top_, epath.c_str()) != -1)
        return;
      sys_fail("Cannot open " + epath);
    }
    sys_fail("Cannot open " + std::string(path));
  }

  //! @brief Construct directory walker.
  //! @param path Directory path or PDB code
  //! @param try_pdbid Character for PDB code expansion
  explicit DirWalk(const std::string& path, char try_pdbid='\0')
    : DirWalk(path.c_str(), try_pdbid) {}

  ~DirWalk() {
    for (auto& d : dirs_)
      tinydir_close(&d.second);
  }

  //! @brief Push directory onto traversal stack.
  //! @param cur_pos Current position in parent directory
  //! @param path Directory path
  //! @throws std::system_error if directory cannot be opened
  void push_dir(size_t cur_pos, const _tinydir_char_t* path) {
    dirs_.emplace_back();
    dirs_.back().first = cur_pos;
    if (tinydir_open_sorted(&dirs_.back().second, path) == -1)
      sys_fail("Cannot open directory " + as_utf8(path));
  }

  //! @brief Pop directory from traversal stack.
  //! @return Position in parent directory
  size_t pop_dir() {
    assert(!dirs_.empty());
    size_t old_pos = dirs_.back().first;
    tinydir_close(&dirs_.back().second);
    dirs_.pop_back();
    return old_pos;
  }

  //! @brief Iterator over files in directory tree.
  struct Iter {
    DirWalk& walk;  //!< Reference to walker
    size_t cur;     //!< Current position in directory

    //! @brief Get current directory.
    //! @return Current directory
    const tinydir_dir& get_dir() const { return walk.dirs_.back().second; }

    //! @brief Get current file.
    //! @return Current file
    const tinydir_file& get() const {
      if (walk.dirs_.empty())
        return walk.top_;
      assert(cur < get_dir().n_files);
      return get_dir()._files[cur];
    }

    //! @brief Dereference to get file path.
    //! @return Current file path
    std::string operator*() const { return as_utf8(get().path); }

    //! @brief Check if name is "." or "..".
    //! @param name File name
    //! @return True if special directory
    bool is_special(const _tinydir_char_t* name) const {
      return name[0] == '.' && (name[1] == '\0' ||
                                (name[1] == '.' && name[2] == '\0'));
    }

    //! @brief Get current directory depth.
    //! @return Depth (0 = top level)
    size_t depth() const { return walk.dirs_.size(); }

    //! @brief Advance to next file (depth-first).
    void next() {
      const tinydir_file& tf = get();
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

    //! @brief Increment to next matching file.
    void operator++() {
      for (;;) {
        next();
        const tinydir_file& f = get();
        if ((!FileOnly && f.is_dir)
            || (!f.is_dir && walk.filter.check(as_utf8(f.name)))
            || walk.is_single_file()
            || (depth() == 0 && cur == 1))
          break;
      }
    }

    //! @brief Compare iterators (used only with end()).
    //! @param o Other iterator
    //! @return True if equal
    bool operator==(const Iter& o) const { return depth()==0 && cur == o.cur; }

    //! @brief Compare iterators.
    //! @param o Other iterator
    //! @return True if not equal
    bool operator!=(const Iter& o) const { return !operator==(o); }
  };

  //! @brief Get begin iterator.
  //! @return Iterator to first matching file
  Iter begin() {
    Iter it{*this, 0};
    if (FileOnly && !is_single_file()) // i.e. the top item is a directory
      ++it;
    return it;
  }

  //! @brief Get end iterator.
  //! @return End iterator
  Iter end() { return Iter{*this, 1}; }

  //! @brief Check if top-level item is a single file.
  //! @return True if single file (not directory)
  bool is_single_file() { return !top_.is_dir; }

private:
  friend struct Iter;
  tinydir_file top_;                               //!< Top-level file/directory
  std::vector<std::pair<size_t, tinydir_dir>> dirs_;  //!< Directory stack
protected:
  Filter filter;  //!< File filter
};

//! Walk over CIF files.
using CifWalk = DirWalk<true, impl::IsCifFile>;

//! Walk over mmCIF files.
using MmCifWalk = DirWalk<true, impl::IsMmCifFile>;

//! Walk over PDB files.
using PdbWalk = DirWalk<true, impl::IsPdbFile>;

//! Walk over coordinate files (excluding SF files).
using CoorFileWalk = DirWalk<true, impl::IsCoordinateFile>;

//! @brief Walk over files matching glob pattern.
struct GlobWalk : public DirWalk<true, impl::IsMatchingFile> {
  //! @brief Construct glob walker.
  //! @param path Directory path
  //! @param glob Glob pattern
  GlobWalk(const std::string& path, const std::string& glob) : DirWalk(path) {
    filter.pattern = glob;
  }
};

} // namespace gemmi
#endif
