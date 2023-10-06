// Copyright 2017 Global Phasing Ltd.

#include <cstdio>
#include <cstdlib>  // for getenv, strtol
#include <cstring>  // for strstr
#include <vector>

#define USE_UNICODE
#ifdef USE_UNICODE
# include <clocale>  // for setlocale
# include <cwchar>   // for wint_t
#endif

inline int terminal_columns() {
  // In bash COLUMNS is a shell variable, not environment variable.
  // It needs to be exported to be visible here.
  if (const char* columns_env = std::getenv("COLUMNS")) {
    long c = std::strtol(columns_env, nullptr, 10);
    if (c > 10 && c < 1000)  // sanity check
      return (int) c;
  }
  return 80;
}

template<typename T>
void print_histogram(const std::vector<T>& data, double min, double max) {
#ifdef USE_UNICODE
  const char* locale = std::setlocale(LC_CTYPE, "");
  const bool use_utf = (locale && std::strstr(locale, "UTF-8") != nullptr);
  const int rows = use_utf ? 12 : 24;
#else
  constexpr int rows = 24;
#endif
  int cols = terminal_columns();
  std::vector<int> bins(cols+1, 0);
  double delta = max - min;
  for (T d : data) {
    int n = (int) std::floor((d - min) * (cols / delta));
    bins[n >= 0 ? (n < cols ? n : cols - 1) : 0]++;
  }
  double max_h = *std::max_element(std::begin(bins), std::end(bins));
  for (int i = rows; i > 0; --i) {
    for (int j = 0; j < cols; ++j) {
      double h = bins[j] / max_h * rows;
#ifdef USE_UNICODE
      if (use_utf) {
        std::wint_t c = ' ';
        if (h > i) {
          c = 0x2588; // 0x2581 = one eighth block, ..., 0x2588 = full block
        } else if (h > i - 1) {
          c = 0x2581 + static_cast<int>((h - (i - 1)) * 7);
        }
        std::printf("%lc", c);
      } else
#endif
      {
        std::putchar(h > i + 0.5 ? '#' : ' ');
      }
    }
    std::putchar('\n');
  }
}

