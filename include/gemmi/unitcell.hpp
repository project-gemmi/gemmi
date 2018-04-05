// Copyright 2017 Global Phasing Ltd.
//
// Unit cell.

#ifndef GEMMI_UNITCELL_HPP_
#define GEMMI_UNITCELL_HPP_

#include <cmath>      // for cos, sin, sqrt, floor
#include <array>
#include "util.hpp"

namespace gemmi {

constexpr double pi() { return 3.1415926535897932384626433832795029; }

struct Vec3 {
  double x, y, z;

  Vec3() : x(0), y(0), z(0) {}
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  double& at(int i) {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default: throw std::out_of_range("Vec3 index must be 0, 1 or 2.");
    }
  }
  double at(int i) const { return const_cast<Vec3*>(this)->at(i); }

  Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
  Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
  Vec3 operator*(double d) const { return {x*d, y*d, z*d}; }
  Vec3 operator/(double d) const { return *this * (1.0/d); }
  Vec3& operator-=(const Vec3& o) { *this = *this - o; return *this; }
  Vec3& operator+=(const Vec3& o) { *this = *this + o; return *this; }
  Vec3& operator*=(double d) { *this = *this * d; return *this; }
  Vec3& operator/=(double d) { return operator*=(1.0/d); }

  Vec3 negated() const { return {-x, -y, -z}; }
  double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
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
  void move_toward_zero_by_one() {
    if (x > 0.5) x -= 1.0; else if (x < -0.5) x += 1.0;
    if (y > 0.5) y -= 1.0; else if (y < -0.5) y += 1.0;
    if (z > 0.5) z -= 1.0; else if (z < -0.5) z += 1.0;
  }
};

enum class SymmetryImage : char { Same, Different, Unspecified };

// Result of find_nearest_image
struct NearbyImage {
  double dist_sq;
  int box[3] = { 0, 0, 0 };
  int sym_id = 0;
  double dist() const { return std::sqrt(dist_sq); }
  bool same_image() const {
    return box[0] == 0 && box[1] == 0 && box[2] == 0 && sym_id == 0;
  }
  std::string pdb_symbol(bool underscore) const {
    char nnn[4] = "555";
    for (int i = 0; i < 3; ++i)
      nnn[i] += box[0];
    return std::to_string(sym_id + 1) + (underscore ? "_" : "") + nnn;
  }
};

struct Matrix33 {
  double a[3][3] = { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} };

  // make it accessible with ".a"
  typedef double row_t[3];
  const row_t& operator[](int i) const { return a[i]; }
  row_t& operator[](int i) { return a[i]; }

  Matrix33() = default;
  Matrix33(double a1, double a2, double a3, double b1, double b2, double b3,
           double c1, double c2, double c3)
  : a{{a1, a2, a3}, {b1, b2, b3}, {c1, c2, c3}} {}

  Vec3 multiply(const Vec3& p) const {
    return {a[0][0] * p.x  + a[0][1] * p.y  + a[0][2] * p.z,
            a[1][0] * p.x  + a[1][1] * p.y  + a[1][2] * p.z,
            a[2][0] * p.x  + a[2][1] * p.y  + a[2][2] * p.z};
  }
  Matrix33 multiply(const Matrix33& b) const {
    Matrix33 r;
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    return r;
  }

  bool approx(const Matrix33& other, double epsilon) const {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if (std::fabs(a[i][j] - other.a[i][j]) > epsilon)
          return false;
    return true;
  }
  double determinant() const {
    return a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]) +
           a[0][1] * (a[1][2]*a[2][0] - a[2][2]*a[1][0]) +
           a[0][2] * (a[1][0]*a[2][1] - a[2][0]*a[1][1]);
  }
  Matrix33 inverse() const {
    Matrix33 inv;
    double inv_det = 1.0 / determinant();
    inv[0][0] = inv_det * (a[1][1] * a[2][2] - a[2][1] * a[1][2]);
    inv[0][1] = inv_det * (a[0][2] * a[2][1] - a[0][1] * a[2][2]);
    inv[0][2] = inv_det * (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
    inv[1][0] = inv_det * (a[1][2] * a[2][0] - a[1][0] * a[2][2]);
    inv[1][1] = inv_det * (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
    inv[1][2] = inv_det * (a[1][0] * a[0][2] - a[0][0] * a[1][2]);
    inv[2][0] = inv_det * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
    inv[2][1] = inv_det * (a[2][0] * a[0][1] - a[0][0] * a[2][1]);
    inv[2][2] = inv_det * (a[0][0] * a[1][1] - a[1][0] * a[0][1]);
    return inv;
  }
  bool is_identity() const {
    return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 &&
           a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 &&
           a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1;
  }
  // Calculate eigenvalues of **symmetric** matrix.
  // Based on https://en.wikipedia.org/wiki/Eigenvalue_algorithm
  std::array<double, 3> calculate_eigenvalues() const {
    double p1 = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
		if (p1 == 0)
      return {{a[0][0], a[1][1], a[2][2]}};
    double q = (1./3.) * (a[0][0] + a[1][1] + a[2][2]);
    Matrix33 b(a[0][0] - q, a[0][1], a[0][2],
               a[1][0], a[1][1] - q, a[1][2],
               a[2][0], a[2][1], a[2][2] - q);
    double p2 = b[0][0] * b[0][0] + b[1][1] * b[1][1] + b[2][2] * b[2][2]
                + 2 * p1;
    double p = std::sqrt((1./6.) * p2);
    double r = b.determinant() / ((1./3.) * p2 * p);
    double phi = 0;
    if (r <= -1)
      phi = (1./3.) * pi();
    else if (r < 1)
      phi = (1./3.) * std::acos(r);
    double eig1 = q + 2 * p * std::cos(phi);
    double eig3 = q + 2 * p * std::cos(phi + 2./3.*pi());
    return {{eig1, 3 * q - eig1 - eig3, eig3}};
  }

  // Assumes symmetric matrix and one of the eigenvalue calculate above.
  // May not work if eigenvalues are not distinct.
  Vec3 calculate_eigenvector(double eigenvalue) const {
    Vec3 r0(a[0][0] - eigenvalue, a[0][1], a[0][2]);
    Vec3 r1(a[1][0], a[1][1] - eigenvalue, a[1][2]);
    Vec3 r2(a[2][0], a[2][1], a[2][2] - eigenvalue);
    Vec3 cr[3] = {r0.cross(r1), r0.cross(r2), r1.cross(r2)};
    int idx = 0;
    double lensq = 0;
    for (int i = 0; i < 3; ++i) {
      double tmp = cr[i].length_sq();
      if (tmp > lensq) {
        idx = i;
        lensq = tmp;
      }
    }
    if (lensq == 0)
      return Vec3(0, 0, 1); // an arbitrary choice for the special case
    return cr[idx] / std::sqrt(lensq);
  }
};

struct Transform {
  Matrix33 mat;
  Vec3 vec;

  Transform inverse() const {
    Matrix33 minv = mat.inverse();
    return {minv, minv.multiply(vec).negated()};
  }

  Vec3 apply(const Vec3& x) const { return mat.multiply(x) + vec; }

  Transform combine(const Transform& b) const {
    return {mat.multiply(b.mat), vec + mat.multiply(b.vec)};
  }

  bool is_identity() const {
    return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
  }
  void set_identity() { mat = Matrix33(); vec = Vec3(); }
};

// for the sake of type safety, a variant that has apply() expecting Fractional
struct FTransform : Transform {
  FTransform(const Transform& t) : Transform(t) {}
  FTransform(Transform&& t) : Transform(t) {}
  FTransform(Matrix33 m, Vec3 v) : Transform{m, v} {}
  Fractional apply(const Fractional& p) const {
    return Fractional(Transform::apply(p));
  }
};


struct UnitCell {
  UnitCell() = default;
  UnitCell(double a_, double b_, double c_,
           double alpha_, double beta_, double gamma_) {
    set(a_, b_, c_, alpha_, beta_, gamma_);
  }
  double a = 1.0, b = 1.0, c = 1.0;
  double alpha = 90.0, beta = 90.0, gamma = 90.0;
  Transform orth;
  Transform frac;
  // volume and reciprocal parameters a*, b*, c*, alpha*, beta*, gamma*
  double volume = 1.0;
  double ar = 1.0, br = 1.0, cr = 1.0;
  double cos_alphar = 0.0, cos_betar = 0.0, cos_gammar = 0.0;
  bool explicit_matrices = false;
  std::vector<FTransform> images;

  // non-crystalline (for example NMR) structures use fake unit cell 1x1x1.
  bool is_crystal() const { return frac.mat[0][0] != 1.0; }

  void calculate_properties() {
    constexpr double deg2rad = pi() / 180.0;
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
    orth.mat = {a,  b * cos_gamma,  c * cos_beta,
                0., b * sin_gamma, -c * cos_alphar_sin_beta,
                0., 0.           ,  c * sin_beta * sin_alphar};
    orth.vec = {0., 0., 0.};

    double o12 = -cos_gamma / (sin_gamma * a);
    double o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma)
                  / (sin_alphar * sin_beta * sin_gamma * a);
    double o23 = cos_alphar / (sin_alphar * sin_gamma * b);
    frac.mat = {1 / a,  o12,                 o13,
                0.,     1 / orth.mat[1][1],  o23,
                0.,     0.,                  1 / orth.mat[2][2]};
    frac.vec = {0., 0., 0.};
  }

  void set_matrices_from_fract(const Transform& f) {
    // mmCIF _atom_sites.fract_transf_* and PDB SCALEn records usually
    // have less significant digits than unit cell parameters, and should
    // be ignored unless we have non-standard settings.
    if (f.mat.approx(frac.mat, 5e-6) && f.vec.approx(frac.vec, 1e-6))
      return;
    // The SCALE record is sometimes incorrect. Here we only catch cases
    // when CRYST1 is set as for non-crystal and SCALE is very suspicious.
    if (frac.mat[0][0] == 1.0 && (f.mat[0][0] == 0.0 || f.mat[0][0] > 1.0))
      return;
    frac = f;
    orth = f.inverse();
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
    return Position(orth.apply(f));
  }
  Fractional fractionalize(const Position& o) const {
    return Fractional(frac.apply(o));
  }

  double volume_per_image() const {
    return is_crystal() ? volume / (1 + images.size()) : NAN;
  }

  // Helper function. PBC = periodic boundary conditions.
  bool search_pbc_images(Fractional&& diff, NearbyImage& image) const {
    int box[3] = { iround(diff.x), iround(diff.y), iround(diff.z) };
    diff.x -= box[0];
    diff.y -= box[1];
    diff.z -= box[2];
    Position orth_diff = orthogonalize(diff);
    double dsq = orth_diff.length_sq();
    if (dsq < image.dist_sq) {
      image.dist_sq = dsq;
      for (int j = 0; j < 3; ++j)
        image.box[j] = box[j];
      return true;
    }
    return false;
  }

  NearbyImage find_nearest_image(const Position& ref, const Position& pos,
                                 SymmetryImage sym_image) const {
    NearbyImage image;
    image.dist_sq = ref.dist_sq(pos);
    if (sym_image == SymmetryImage::Same || !is_crystal()) {
      if (sym_image == SymmetryImage::Different || image.dist_sq == 0.0)
        image.dist_sq = INFINITY;
      return image;
    }
    Fractional fpos = fractionalize(pos);
    Fractional fref = fractionalize(ref);
    search_pbc_images(fpos - fref, image);
    if ((sym_image == SymmetryImage::Different || image.dist_sq == 0.0) &&
        image.box[0] == 0 && image.box[1] == 0 && image.box[2] == 0)
      image.dist_sq = INFINITY;
    for (int n = 0; n != static_cast<int>(images.size()); ++n)
      if (search_pbc_images(images[n].apply(fpos) - fref, image))
        image.sym_id = n + 1;
    return image;
  }

  // return number of nearby symmetry mates (0 = none, 3 = 4-fold axis, etc)
  int is_special_position(const Position& pos, double max_dist = 0.8) const {
    const double max_dist_sq = max_dist * max_dist;
    int n = 0;
    Fractional fpos = fractionalize(pos);
    for (const FTransform& image : images)
      if (orthogonalize(image.apply(fpos) - fpos).length_sq() < max_dist_sq)
        ++n;
    return n;
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
