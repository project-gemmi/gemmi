// Copyright Global Phasing Ltd.
//
// Logger - a tiny utility for passing messages through a callback.

#ifndef GEMMI_LOGGER_HPP_
#define GEMMI_LOGGER_HPP_

#include <cstdio>      // for fprintf
#include <functional>  // for function
#include "fail.hpp"    // for GEMMI_COLD
#include "util.hpp"    // for cat

namespace gemmi {

/// Passes notes and warnings to a callback function.
/// Note: if the callback is not set, notes are ignored but warnings
/// are thrown as exceptions.
struct Logger {
  using Callback = std::function<void(const std::string&)>;
  Callback callback;

  // For internal use in functions that produce messages: suspending when
  // the same function is called multiple times avoids duplicated messages.
  bool suspended = false;

  // Send warning.
  template<class... Args> GEMMI_COLD void err(Args const&... args) const {
    if (!suspended) {
      std::string msg = cat(args...);
      if (callback == nullptr)
        fail(msg);
      callback("Warning: " + msg);
    }
  }

  // Send note.
  template<class... Args> void note(Args const&... args) const {
    if (!suspended && callback)
      callback(cat("Note: ", args...));
  }

  // Send a message without any prefix
  template<class... Args> void mesg(Args const&... args) const {
    if (!suspended && callback)
      callback(cat(args...));
  }

  /// to be used as: logger.callback = Logger::to_stderr;
  static void to_stderr(const std::string& s) {
    std::fprintf(stderr, "%s\n", s.c_str());
  }
  /// to be used as: logger.callback = Logger::to_stdout;
  static void to_stdout(const std::string& s) {
    std::fprintf(stdout, "%s\n", s.c_str());
  }
  /// to be used as: logger.callback = Logger::nop;
  static void nop(const std::string&) {}
};

} // namespace gemmi
#endif
