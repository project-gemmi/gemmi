// Copyright 2017 Global Phasing Ltd.
//
// Unit cell.

#ifndef GEMMI_UNITCELL_HPP_
#define GEMMI_UNITCELL_HPP_

#include <cmath>      // for cos, sin, sqrt, floor
#include <linalg.h>
#include "util.hpp"

namespace gemmi {

struct Vec3 {
  //Vec3() {}
  //Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  double x, y, z;
  // FIXME: it may be UB, switch to array and references x, y, z
  constexpr double operator[](int i) const { return (&x)[i]; }
  double& operator[](int i) { return (&x)[i]; }
  Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
  Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
  double dot(const Vec3& o) const { return x*o.x + y*o.x + z*o.z; }
  Vec3 cross(const Vec3& o) const {
    return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
  }
  double length_sq() const { return x * x + y * y + z * z; }
  double dist_sq(const Vec3& o) const { return (*this - o).length_sq(); }
  double dist(const Vec3& o) const { return std::sqrt(dist_sq(o)); }
  bool approx(const Vec3& o, double epsilon) const {
    return std::fabs(x - o.x) <= epsilon &&
           std::fabs(y - o.y) <= epsilon &&
           std::fabs(z - o.z) <= epsilon;
  }
  std::string str() const {
    using namespace std;
    char buf[64] = {0};
    snprintf(buf, 63, "[%g %g %g]", x, y, z);
    return buf;
  }
};

// coordinates in Angstroms (a.k.a. orthogonal coordinates)
struct Position : Vec3 {
  Position() = default;
  Position(double x_, double y_, double z_) : Vec3{x_, y_, z_} {}
  explicit Position(Vec3&& v) : Vec3(v) {}
  Position operator-(const Position& o) const {
    return Position(Vec3::operator-(o));
  }
  Position operator+(const Position& o) const {
    return Position(Vec3::operator+(o));
  }
};

// fractional coordinates
struct Fractional : Vec3 {
  Fractional() = default;
  Fractional(double x_, double y_, double z_) : Vec3{x_, y_, z_} {}
  explicit Fractional(Vec3&& v) : Vec3(v) {}
  Fractional operator-(const Fractional& o) const {
    return Fractional(Vec3::operator-(o));
  }
  Fractional operator+(const Fractional& o) const {
    return Fractional(Vec3::operator+(o));
  }
  Fractional& wrap_to_unit() {
    x -= std::floor(x);
    y -= std::floor(y);
    z -= std::floor(z);
    return *this;
  }
};

// Result of find_nearest_image
struct NearestImage {
  double dist_sq;
  int box[3] = { 0, 0, 0 };
  int sym_id = 0;
  std::string pdb_symbol(bool underscore) const {
    char nnn[4] = "555";
    for (int i = 0; i < 3; ++i)
      nnn[i] += box[0];
    return std::to_string(sym_id + 1) + (underscore ? "_" : "") + nnn;
  }
};

struct Matrix33 {
  double a[3][3] = { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} };

  Matrix33() = default;
  Matrix33(double a1, double a2, double a3, double b1, double b2, double b3,
           double c1, double c2, double c3)
  : a{{a1, a2, a3}, {b1, b2, b3}, {c1, c2, c3}} {}

  Vec3 multiply(const Vec3& p) const {
    return {a[0][0] * p.x  + a[0][1] * p.y  + a[0][2] * p.z,
            a[1][0] * p.x  + a[1][1] * p.y  + a[1][2] * p.z,
            a[2][0] * p.x  + a[2][1] * p.y  + a[2][2] * p.z};
  }

  bool approx(const Matrix33& other, double epsilon) const {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if (std::fabs(a[i][j] - other.a[i][j]) > epsilon)
          return false;
    return true;
  }
};

struct Matrix44 {
  double a[4][4] = {{1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}};
  bool is_identity() const {
    return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 && a[0][3] == 0 &&
           a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 && a[1][3] == 0 &&
           a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1 && a[2][3] == 0 &&
           a[3][0] == 0 && a[3][1] == 0 && a[3][2] == 0 && a[3][3] == 1;
  }
  Matrix33 mat33() const {
    return Matrix33(a[0][0], a[0][1], a[0][2],
                    a[1][0], a[1][1], a[1][2],
                    a[2][0], a[2][1], a[2][2]);
  }
  Vec3 shift() const { return { a[3][0], a[3][1], a[3][2] }; }
  double determinant() const {
      return a[0][0]*(a[1][1]*a[2][2]*a[3][3] + a[3][1]*a[1][2]*a[2][3] +
                      a[2][1]*a[3][2]*a[1][3] - a[1][1]*a[3][2]*a[2][3] -
                      a[2][1]*a[1][2]*a[3][3] - a[3][1]*a[2][2]*a[1][3]) +
             a[0][1]*(a[1][2]*a[3][3]*a[2][0] + a[2][2]*a[1][3]*a[3][0] +
                      a[3][2]*a[2][3]*a[1][0] - a[1][2]*a[2][3]*a[3][0] -
                      a[3][2]*a[1][3]*a[2][0] - a[2][2]*a[3][3]*a[1][0]) +
             a[0][2]*(a[1][3]*a[2][0]*a[3][1] + a[3][3]*a[1][0]*a[2][1] +
                      a[2][3]*a[3][0]*a[1][1] - a[1][3]*a[3][0]*a[2][1] -
                      a[2][3]*a[1][0]*a[3][1] - a[3][3]*a[2][0]*a[1][1]) +
             a[0][3]*(a[1][0]*a[3][1]*a[2][2] + a[2][0]*a[1][1]*a[3][2] +
                      a[3][0]*a[2][1]*a[1][2] - a[1][0]*a[2][1]*a[3][2] -
                      a[3][0]*a[1][1]*a[2][2] - a[2][0]*a[3][1]*a[1][2]);
  }
  Matrix44 inverse() const {
    // TODO
    return *this;
  }
};

// discussion: https://stackoverflow.com/questions/20305272/
inline double calculate_dihedral(const Position& p0, const Position& p1,
                                 const Position& p2, const Position& p3) {
  Vec3 b0 = p1 - p0;
  Vec3 b1 = p2 - p1;
  Vec3 b2 = p3 - p2;
  Vec3 u = b1.cross(b0);
  Vec3 w = b2.cross(b1);
  double y = u.cross(w).dot(b1);
  double x = u.dot(w) * std::sqrt(b1.length_sq());
  return std::atan2(y, x);
}


struct UnitCell {
  double a = 1.0, b = 1.0, c = 1.0;
  double alpha = 90.0, beta = 90.0, gamma = 90.0;
  constexpr double operator[](int i) const { return (&a)[i]; }
  Matrix33 orth;
  Matrix33 frac;
  Vec3 shift = {0., 0., 0.};
  // volume and reciprocal parameters a*, b*, c*, alpha*, beta*, gamma*
  double volume = 1.0;
  double ar = 1.0, br = 1.0, cr = 1.0;
  double cos_alphar = 0.0, cos_betar = 0.0, cos_gammar = 0.0;
  bool explicit_matrices = false;

  // non-crystalline (for example NMR) structures use fake unit cell 1x1x1.
  bool is_crystal() const { return a != 1.0; }

  void calculate_properties() {
    constexpr double deg2rad = 3.1415926535897932384626433832795029 / 180.0;
    // ensure exact values for right angles
    double cos_alpha = alpha == 90. ? 0. : std::cos(deg2rad * alpha);
    double cos_beta  = beta  == 90. ? 0. : std::cos(deg2rad * beta);
    double cos_gamma = gamma == 90. ? 0. : std::cos(deg2rad * gamma);
    double sin_alpha = alpha == 90. ? 1. : std::sin(deg2rad * alpha);
    double sin_beta  = beta  == 90. ? 1. : std::sin(deg2rad * beta);
    double sin_gamma = gamma == 90. ? 1. : std::sin(deg2rad * gamma);
    if (sin_alpha == 0 || sin_beta == 0 || sin_gamma == 0)
      gemmi::fail("Impossible angle - N*180deg.");

		// volume - formula from Giacovazzo p.62
    volume = a * b * c * sqrt(1 - cos_alpha * cos_alpha - cos_beta * cos_beta
                              - cos_gamma * cos_gamma
                              + 2 * cos_alpha * cos_beta * cos_gamma);

    // reciprocal parameters a*, b*, ... (Giacovazzo, p. 64)
    ar = b * c * sin_alpha / volume;
    br = a * c * sin_beta / volume;
    cr = a * b * sin_gamma / volume;
    double cos_alphar_sin_beta = (cos_beta * cos_gamma - cos_alpha) / sin_gamma;
    cos_alphar = cos_alphar_sin_beta / sin_beta;
    //cos_alphar = (cos_beta * cos_gamma - cos_alpha) / (sin_beta * sin_gamma);
    cos_betar = (cos_alpha * cos_gamma - cos_beta) / (sin_alpha * sin_gamma);
    cos_gammar = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta);

    if (explicit_matrices)
      return;

    // The orthogonalization matrix we use is described in ITfC B p.262:
    // "An alternative mode of orthogonalization, used by the Protein
    // Data Bank and most programs, is to align the a1 axis of the unit
    // cell with the Cartesian X_1 axis, and to align the a*_3 axis with the
    // Cartesian X_3 axis."
    double sin_alphar = std::sqrt(1.0 - cos_alphar * cos_alphar);
    orth = {a,  b * cos_gamma,  c * cos_beta,
            0., b * sin_gamma, -c * cos_alphar_sin_beta,
            0., 0.           ,  c * sin_beta * sin_alphar};

    double o12 = -cos_gamma / (sin_gamma * a);
    double o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma)
                  / (sin_alphar * sin_beta * sin_gamma * a);
    double o23 = cos_alphar / (sin_alphar * sin_gamma * b);
    frac = {1 / a,  o12,               o13,
            0.,     1 / orth.a[1][1],  o23,
            0.,     0.,                1 / orth.a[2][2]};
  }

  void set_matrices_from_fract(const linalg::mat<double,4,4>& fract) {
    // mmCIF _atom_sites.fract_transf_* and PDB SCALEn records usually
    // have less significant digits than unit cell parameters, and should
    // be ignored unless we have non-standard settings.
    Matrix33 f;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        f.a[i][j] = fract[j][i];
    Vec3 s = {fract.w.x, fract.w.y, fract.w.z};
    if (f.approx(frac, 1e-6) && s.approx(shift, 1e-6))
      return;
    frac = f;
    shift = s;
    auto ortho = linalg::inverse(fract);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        orth.a[i][j] = ortho[j][i];
    explicit_matrices = true;
  }

  void set_matrices_from_fract(const Matrix44& f) {
    frac = f.mat33();
    shift = f.shift();
    orth = f.inverse().mat33();
    explicit_matrices = true;
  }

  void set(double a_, double b_, double c_,
           double alpha_, double beta_, double gamma_) {
    if (gamma_ == 0.0)  // ignore empty/partial CRYST1 (example: 3iyp)
      return;
    a = a_;
    b = b_;
    c = c_;
    alpha = alpha_;
    beta = beta_;
    gamma = gamma_;
    calculate_properties();
  }

  // we could also apply shift for the few special cases that have
  // SCALEn with non-zero vector
  Position orthogonalize(const Fractional& f) const {
    return Position(orth.multiply(f));
  }
  Fractional fractionalize(const Position& o) const {
    return Fractional(frac.multiply(o));
  }
};


struct SymmetryOp {
  Matrix33 rot;
  Vec3 tran;
  Fractional apply(const Fractional& p) const {
    return Fractional(rot.multiply(p) + tran);
  }
};

struct UnitCellWithSymmetry : UnitCell {
  std::vector<SymmetryOp> images;

  // Helper function. PBC = periodic boundary conditions.
  void search_pbc_images(Fractional&& diff, NearestImage& image, int id) const {
    int box[3];
    for (int i = 0; i < 3; ++i) {
      box[i] = iround(diff[i]);
      diff[i] -= box[i];
    }
    Position orth_diff = orthogonalize(diff);
    double dsq = orth_diff.length_sq();
    if (dsq < image.dist_sq) {
      image.dist_sq = dsq;
      for (int j = 0; j < 3; ++j)
        image.box[j] = box[j];
      image.sym_id = id;
    }
    //std::fprintf(stderr, " [%d] d = %g  %d %d %d\n", id, std::sqrt(dsq),
    //                     box[0], box[1], box[2]);
  }

  NearestImage find_nearest_image(const Position& ref, const Position& pos,
                                  bool non_ident) const {
    NearestImage image;
    image.dist_sq = ref.dist_sq(pos);
    if (!is_crystal())
      return image;
    Fractional fpos = fractionalize(pos);
    Fractional fref = fractionalize(ref);
    search_pbc_images(fpos - fref, image, 0);
    if (non_ident && image.dist_sq == 0 &&
        image.box[0] == 0 && image.box[1] == 0 && image.box[2] == 0)
      image.dist_sq = 1e100;
    for (int n = 0; n != static_cast<int>(images.size()); ++n)
      search_pbc_images(images[n].apply(fpos) - fref, image, n+1);
    return image;
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
