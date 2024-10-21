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
/// Messages are passed as strings without a newline character.
/// Messages have severity levels similar syslog:
///  7=debug, 6=info (all but debug), 5=notice, 3=error
/// A numeric threshold can be set to limit the messages (see below).
/// Quirk: if a callback is not set, errors are thrown as exceptions.
struct Logger {
  /// A function that handles messages.
  std::function<void(const std::string&)> callback;
  /// Pass messages of this level and all lower (more severe) levels:
  /// 7=all messages, 6=all but debug, 0=none
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

  /// Send a warning/error. Unrecoverable errors are thrown directly and
  /// don't go through this class, so here we're left with errors that
  /// can be downgraded to warnings. If a callback is set, the message is
  /// passed as a warning; otherwise it's thrown as a std::runtime_error.
  /// (Admittedly, it's a questionable design.)
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
  /// to be used as: logger.callback = Logger::nop;
  static void nop(const std::string&) {}
};

} // namespace gemmi
#endif
