// Code in this file is based on and derived from files ksw2_gg.c and ksw2.h
// from https://github.com/lh3/ksw2, which is under the MIT license.
// The original code, written by Heng Li, has more features and has more
// efficient variants that use SSE instructions.

#ifndef GEMMI_SEQALIGN_HPP_
#define GEMMI_SEQALIGN_HPP_

#include <cstdint>
#include <algorithm> // for reverse
#include <string>
#include <vector>

namespace gemmi {

struct Alignment {
  struct Item {
    std::uint32_t value;
    char op() const { return "MID"[value & 0xf]; }
    std::uint32_t len() const { return value >> 4; }
  };
  int score;
  std::vector<Item> cigar;

  std::string cigar_str() const {
    std::string s;
    for (Item item : cigar) {
      s += std::to_string(item.len());
      s += item.op();
    }
    return s;
  }

  // In the backtrack matrix, value p[] has the following structure:
  //   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F
  //   bit 3/0x08: 1 if a continuation on the E state
  //   bit 4/0x10: 1 if a continuation on the F state
  void backtrack_to_cigar(const std::uint8_t *p, int i, int j) {
    i--;
    int j0 = j--;
    int state = 0;
    while (i >= 0 && j >= 0) {
      // at the beginning of the loop, _state_ tells us which state to check
      // if requesting the H state, find state one maximizes it.
      uint32_t tmp = p[(size_t)i * j0 + j];
      if (state == 0 || (tmp & (1 << (state + 2))) == 0)
        state = tmp & 7;
      if (state == 0) { // match
        push_cigar(0, 1);
        --i;
        --j;
      } else if (state == 1) { // deletion
        push_cigar(2, 1);
        --i;
      } else { // insertion
        push_cigar(1, 1);
        --j;
      }
    }
    if (i >= 0)
      push_cigar(2, i + 1); // first deletion
    else if (j >= 0)
      push_cigar(1, j + 1); // first insertion
    std::reverse(cigar.begin(), cigar.end());
  }

private:
  void push_cigar(std::uint32_t op, int len) {
    if (cigar.empty() || op != (cigar.back().value & 0xf))
      cigar.push_back({len<<4 | op});
    else
      cigar.back().value += len<<4;
  }
};

inline
Alignment align_sequences(int qlen, const std::uint8_t *query,
                          int tlen, const std::uint8_t *target,
                          const std::vector<bool>& free_gapo,
                          std::int8_t m, const std::int8_t *mat,
                          std::int8_t gapo, std::int8_t gape) {
  // generate the query profile
  std::int8_t *query_profile = new std::int8_t[qlen * m];
  for (std::int32_t k = 0, i = 0; k < m; ++k)
    for (std::int32_t j = 0; j < qlen; ++j)
      query_profile[i++] = mat[k * m + query[j]];

  struct eh_t { std::int32_t h, e; };
  eh_t *eh = new eh_t[qlen + 1];
  std::int32_t gapoe = gapo + gape;

  // fill the first row
  {
    std::int32_t gap0 = !free_gapo.empty() && free_gapo[0] ? gape : gapoe;
    eh[0].h = 0;
    eh[0].e = -gap0 - gapoe;
    for (std::int32_t j = 1; j <= qlen; ++j) {
      eh[j].h = -(gap0 + gape * (j - 1));
      eh[j].e = -(gap0 + gapoe + gape * j);
    }
  }

  // backtrack matrix; in each cell: f<<4|e<<2|h
  std::uint8_t *z = new std::uint8_t[(size_t)qlen * tlen];
  // DP loop
  for (std::int32_t i = 0; i < tlen; ++i) {
    std::uint8_t target_item = target[i];
    std::int8_t *scores = &query_profile[target_item * qlen];
    std::uint8_t *zi = &z[(size_t)i * qlen];
    std::int32_t h1 = -(gapoe + gape * i);
    std::int32_t f = -(gapoe + gapoe + gape * i);
    std::int32_t gapx = i+1 < (std::int32_t)free_gapo.size() && free_gapo[i+1]
                        ? gape : gapoe;
    for (std::int32_t j = 0; j < qlen; ++j) {
      // At the beginning of the loop:
      //  eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
      // Cells are computed in the following order:
      //   H(i,j)   = max{H(i-1,j-1) + S(i,j), E(i,j), F(i,j)}
      //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
      //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
      eh_t *p = &eh[j];
      std::int32_t h = p->h;
      std::int32_t e = p->e;
      p->h = h1;
      h += scores[j];
      std::uint8_t direction = 0;
      if (h < e) {
        direction = 1;  // deletion
        h = e;
      }
      if (h <= f) {
        direction = 2;  // insertion
        h = f;
      }
      h1 = h;

      h -= gapoe;
      e -= gape;
      if (e > h)
        direction |= 0x08;
      else
        e = h;

      h = h1 - gapx;
      p->e = e;
      f -= gape;
      if (f > h)
        direction |= 0x10;
      else
        f = h;

      // z[i,j] keeps h for the current cell and e/f for the next cell
      zi[j] = direction;
    }
    eh[qlen].h = h1;
    eh[qlen].e = -0x40000000; // -infinity
  }

  Alignment result;
  result.score = eh[qlen].h;
  delete [] query_profile;
  delete [] eh;
  result.backtrack_to_cigar(z, tlen, qlen);
  delete [] z;
  return result;
}

} // namespace gemmi
#endif
