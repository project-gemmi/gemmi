// Copyright 2018 Global Phasing Ltd.
//
// Math utilities. 3D linear algebra.

#ifndef GEMMI_MATH_HPP_
#define GEMMI_MATH_HPP_

#include <cmath>      // for fabs, cos, sqrt, round
#include <array>
#include <stdexcept>  // for out_of_range
#include <string>
#include <cstdio>     // for snprintf

namespace gemmi {

constexpr double pi() { return 3.1415926535897932384626433832795029; }

inline double deg(double angle) { return 180.0 / pi() * angle; }

inline float sq(float x) { return x * x; }
inline double sq(double x) { return x * x; }

inline int iround(double d) { return static_cast<int>(std::round(d)); }

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
  double length() const { return std::sqrt(length_sq()); }
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

struct Mat33 {
  double a[3][3] = { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} };

  // make it accessible with ".a"
  typedef double row_t[3];
  const row_t& operator[](int i) const { return a[i]; }
  row_t& operator[](int i) { return a[i]; }

  Mat33() = default;
  Mat33(double a1, double a2, double a3, double b1, double b2, double b3,
        double c1, double c2, double c3)
  : a{{a1, a2, a3}, {b1, b2, b3}, {c1, c2, c3}} {}

  Vec3 multiply(const Vec3& p) const {
    return {a[0][0] * p.x  + a[0][1] * p.y  + a[0][2] * p.z,
            a[1][0] * p.x  + a[1][1] * p.y  + a[1][2] * p.z,
            a[2][0] * p.x  + a[2][1] * p.y  + a[2][2] * p.z};
  }
  Mat33 multiply(const Mat33& b) const {
    Mat33 r;
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    return r;
  }
  Mat33 transpose() const {
    return Mat33(a[0][0], a[1][0], a[2][0],
                 a[0][1], a[1][1], a[2][1],
                 a[0][2], a[1][2], a[2][2]);
  }

  bool approx(const Mat33& other, double epsilon) const {
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
  Mat33 inverse() const {
    Mat33 inv;
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
    Mat33 b(a[0][0] - q, a[0][1], a[0][2],
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
  Mat33 mat;
  Vec3 vec;

  Transform inverse() const {
    Mat33 minv = mat.inverse();
    return {minv, minv.multiply(vec).negated()};
  }

  Vec3 apply(const Vec3& x) const { return mat.multiply(x) + vec; }

  Transform combine(const Transform& b) const {
    return {mat.multiply(b.mat), vec + mat.multiply(b.vec)};
  }

  bool is_identity() const {
    return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
  }
  void set_identity() { mat = Mat33(); vec = Vec3(); }
};

struct BoundingBox {
  Vec3 low = Vec3(INFINITY, INFINITY, INFINITY);
  Vec3 high = Vec3(-INFINITY, -INFINITY, -INFINITY);
  void add(const Vec3& p) {
    if (p.x < low.x) low.x = p.x; if (p.x > high.x) high.x = p.x;
    if (p.y < low.y) low.y = p.y; if (p.y > high.y) high.y = p.y;
    if (p.z < low.z) low.z = p.z; if (p.z > high.z) high.z = p.z;
  }
  Vec3 get_size() const { return high - low; }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
