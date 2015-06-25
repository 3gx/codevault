#pragma once
#include <stdexcept>

template<typename OS>
static void safe_printf(OS &os, const char *s)
{
  while (*s) {
    if (*s == '%') {
      if (*(s + 1) == '%') {
        ++s;
      }
      else {
        throw std::runtime_error("invalid format string: missing arguments");
      }
    }
    os << *s++;
  }
}

template<typename OS, typename T, typename... Args>
static void safe_printf(OS &os, const char *s, T&& value, Args&&... args)
{
  while (*s) {
    if (*s == '%') {
      if (*(s + 1) == '%') {
        ++s;
      }
      else {
        os << value;
        safe_printf(os, s + 1, args...); // call even when *s == 0 to detect extra arguments
        return;
      }
    }
    os << *s++;
  }
  throw std::runtime_error("extra arguments provided to safe_printf");
}

