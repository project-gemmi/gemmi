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

/// Passes messages (including warnings/errors) to a callback function.
/// Messages are passed as strings without a trailing newline.
/// They have syslog-like severity levels: 7=debug, 6=info, 5=notice, 3=error,
/// allowing the use of a threshold to filter them.
/// Quirk: Errors double as both errors and warnings. Unrecoverable errors
///        don't go through this class; Logger only handles errors that can
///        be downgraded to warnings. If a callback is set, the error is passed
///        as a warning message. Otherwise, it's thrown as std::runtime_error.
struct Logger {
  /// A function that handles messages.
  std::function<void(const std::string&)> callback;
  /// Pass messages of this level and all lower (more severe) levels:
  /// 7=all, 6=all but debug, 5=notes and warnings, 3=warnings, 0=none
  int threshold = 6;

  /// suspend() and resume() are used internally to avoid duplicate messages
  /// when the same function is called (internally) multiple times.
  void suspend() { threshold -= 100; }
  void resume()  { threshold += 100; }

  /// Send a debug message.
  template<class... Args> void debug(Args const&... args) const {
    if (threshold >= 7 && callback)
      callback(cat("Debug: ", args...));
  }

  /// Send a message without any prefix.
  template<class... Args> void mesg(Args const&... args) const {
    if (threshold >= 6 && callback)
      callback(cat(args...));
  }

  /// Send a note (a notice, a significant message).
  template<class... Args> void note(Args const&... args) const {
    if (threshold >= 5 && callback)
      callback(cat("Note: ", args...));
  }

  /// Send a warning/error (see Quirk above).
  template<class... Args> GEMMI_COLD void err(Args const&... args) const {
    if (threshold >= 3) {
      std::string msg = cat(args...);
      if (callback == nullptr)
        fail(msg);
      callback("Warning: " + msg);
    }
  }

  // predefined callbacks

  /// to be used as: logger.callback = Logger::to_stderr;
  static void to_stderr(const std::string& s) {
    std::fprintf(stderr, "%s\n", s.c_str());
  }
  /// to be used as: logger.callback = Logger::to_stdout;
  static void to_stdout(const std::string& s) {
    std::fprintf(stdout, "%s\n", s.c_str());
  }
};

} // namespace gemmi
#endif
