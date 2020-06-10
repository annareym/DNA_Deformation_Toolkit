#ifndef DnaParams_h
#define DnaParams_h

#include <Eigen/Core>
#include <Eigen/Geometry>  // cross-product
#include <cmath>
#include <cassert>  // assert()
#include <ostream>
#include <map>

namespace DnaParams {
namespace Detail {

/*
DNA parameters according to the scheme described by Lavery et al. [LMMPZ09].

Alternative scheme is by Olson et al. [OBBD+01].

References:
    
    [LMMPZ09]. Lavery R, Moakher M, Maddocks JH, Petkeviciute D, Zakrzewska K.
    Conformational analysis of nucleic acids revisited: Curves+. 
    Nucleic Acids Res. 2009;37(17):5917‐5929. doi:10.1093/nar/gkp608
    
    [OBBD+01]. Olson WK, Bansal M, Burley SK, Dickerson RE, et al.
    A standard reference frame for the description of nucleic acid base-pair geometry.
    J Mol Biol. 2001;313(1):229‐237. doi:10.1006/jmbi.2001.4987
*/

const double PI = 3.141592;

// Tsukuba DNA Base Geometry Parameters for "[LMMPZ09]"-frame.
const double d = 4.702;
const double tau_1 = 2.4691;
const double tau_2 = -0.95138;

template <typename T>
using Matrix = Eigen::Matrix<T,3,3>;

template <typename T>
using Vector = Eigen::Matrix<T,3,1>;


template <typename T>
struct Frame {
    Matrix<T> basis;
    Vector<T> origin;
};

template <typename T>
struct AxisAngle {
    Vector<T> axis;
    T angle;
};

template <typename T> 
static inline void assert_unitary(const Matrix<T> &M) {
    //assert(M.isUnitary());
}

template <typename T>
std::ostream &operator<<(std::ostream &os, Frame<T> const &f) { 
    return os << "basis:" << std::endl << f.basis << std::endl << "origin:" 
    << std::endl << f.origin << std::endl;
}


// Returns vector obtained by rotating vector `v` around `axis` on `angle` radians.
template <typename T>
static inline Vector<T> rotateRodrigues(const Vector<T> &axis, const T angle, const Vector<T> &v) {
    using std::cos;
    return v * cos(angle) + axis.cross(v) * sin(angle) + \
        axis*(axis.dot(v))*(1-cos(angle));
}

template <typename T>
static inline Vector<T> rotateRodrigues(const AxisAngle<T> &axisAngle, const Vector<T> &v) {
    return rotateRodrigues(axisAngle.axis, axisAngle.angle, v);
}

// Returns coordinate frame (triad) using algorithm from [LMMPZ09].
// 
// [LMMPZ09]: In order to avoid having to give the reference system in
// Cartesian coordinates for each standard base, we calculate it using
// chosen base atoms. These are C1', N1(Y)/N9(R) and C2(Y)/C4(R) in standard
// bases (where Y is a pyrimidine and R is a purine).
template <typename T>
static inline Frame<T> frameSb(const Vector<T> &C1p, const Vector<T> &N, const Vector<T> &C) {
    // [LMMPZ09]: The direction of the normal b_N is given by the cross
    // product (N1-C1') x (N1-C2) for pyrimidines and (N9-C1') x (N9-C4) for
    // purines.
    const Vector<T> b_N(((N - C1p).cross(N - C)).normalized());  // "Direction Z".

    // [LMMPZ09]: The base reference point (termed b_R below) is obtained by
    // rotating a vector of length $d$ (initially aligned with the N-C1'
    // direction) clockwise by an angle tau1 around the normal vector passing
    // through the N atom.
    // Comment:  "N-C1' direction" appears to mean from N to C1', i.e. dir=C1'-N.
    const Vector<T> dir_n_c1p((C1p - N).normalized());  // Used for b_R and b_L.
    const Vector<T> v(T(d) * dir_n_c1p);
    const Vector<T> b_R(rotateRodrigues<T>(b_N, tau_1, v) + N);  // "Origin".
    
    // [LMMPZ09]: The next vector of the reference system, pointing towards
    // the phosphodiester backbone joined to the base (termed b_L below) is
    // obtained by a similar rotation, but using a unit vector and the angle
    // tau2.
    const Vector<T> b_L(rotateRodrigues<T>(b_N, tau_2, dir_n_c1p));  // "Direction Y".

    // [LMMPZ09]: The last vector of the reference system, pointing into the
    // major groove, b_D, is obtained from the cross product $b_L x b_N$.
    const Vector<T> b_D(b_L.cross(b_N));  // "Direction X".

    Matrix<T> basis; basis << b_D, b_L, b_N;  // x, y, z,
    assert_unitary(basis);

    const Frame<T> frame = { .basis = basis, .origin = b_R };
    return frame;
}

// Returns rotation angle for a rotation matrix.
// Formula: trace(Q) = 1 + 2 * cos(angle)
// Rotation matrix have specific egenvalues: 1, exp(i*angle), exp(-i*angle).
// Trace of a matrix equals the sum of its eigenvalues.
template <typename T>
static inline T angleFromMatrix(const Matrix<T> &Q) {
    assert_unitary(Q);
    using std::acos;
    const T angle = acos((Q.trace() - 1.0)/2.0);
    return angle;
}

// Returns rotation axis from a rotation matrix.
// See e.g. https://stackoverflow.com/a/12472591
template <typename T>
static inline Vector<T> axisFromMatrix(const Matrix<T> &Q) {
    assert_unitary(Q);
    const T x = Q(2,1) - Q(1,2);
    const T y = Q(0,2) - Q(2,0);
    const T z = Q(1,0) - Q(0,1);
    Vector<T> axis; axis << x, y, z;
    axis.normalize();
    return axis;
}

template <typename T>
static inline AxisAngle<T> axisAngleFromMatrix(const Matrix<T> &Q) {
    return { .axis  = axisFromMatrix(Q),
             .angle = angleFromMatrix(Q) };
}


// Returns rotation matrix that transforms frame cf1 into cf2.
// 
// Frames cf1 and cf2 consist of three column-vectors stacked horizontally:
//       cf1 = [a b c], cf2 = [d e f].
// Algorithm:
//       Q * [a b c] = [d e f]  =>  Q = [d e f] * inv([a b c])
// Instead of inversion, we can use transpose, since [a b c] are orthogonal:
//       inv([a b c]) = [a b c]'.
template <typename T>
static inline Matrix<T> rotmatFromFrames(const Frame<T> &cf1, const Frame<T> &cf2) {
    const Matrix<T> abc(cf1.basis);
    const Matrix<T> def(cf2.basis);
    assert_unitary(abc);
    assert_unitary(def);    
    const Matrix<T> Q(def * abc.transpose());
    assert_unitary(Q);

    return Q;
}


template <typename T>
static inline Matrix<T> rotated(const AxisAngle<T> axisAngle, const Matrix<T> &M) {
    const Vector<T> x(M.col(0));
    const Vector<T> y(M.col(1));
    const Vector<T> z(M.col(2));

    const Vector<T> xm(rotateRodrigues(axisAngle, x));
    const Vector<T> ym(rotateRodrigues(axisAngle, y));
    const Vector<T> zm(rotateRodrigues(axisAngle, z));

    Matrix<T> result; result << xm, ym, zm;    
    return result;
}

template <typename T>
static inline Matrix<T> rotatedHalfway(const Vector<T> &axis, const T &angle, const Matrix<T> &M) {
    const AxisAngle<T> halfRotation = { 
        .axis = axis,
        .angle = angle / 2.0,
    };
    return rotated(halfRotation, M);
}

// Rotate basis A half-way to B.
template <typename T>
static inline Matrix<T> midbasis(const Frame<T> &A, const Frame<T> &B) {
    const Matrix<T> Q(rotmatFromFrames(A, B));
    const auto [axis, angle] = axisAngleFromMatrix(Q);

    const Matrix<T> Mb(rotatedHalfway(axis, angle, A.basis));
    assert_unitary(Mb);
    return Mb;
}

// Returns a "mid-frame": half-way rotated frame from frame A to B.
//
// [LMMPZ09]: It is convenient to express [base-pair and step parameters] 
// with respect to components in a mean reference. To do this as symmetrically 
// as possible, we choose an average frame that is obtained by rotation and 
// translation of the first base reference system, but now through the half
// angle, about the same axis vector, and with the half translation.
template <typename T>
static inline Frame<T> midframe(const Frame<T> &A, const Frame<T> &B) {
    return { 
        .basis  = midbasis(A, B), 
        .origin = (A.origin + B.origin) / T(2.0) 
    };
}

// Returns base-pair coordinate frame from two DNA-base coordinate frame.
template <typename T>
static inline Frame<T> frameBp(const Frame<T> &A, const Frame<T> &B) {
    // [LMMPZ09]: Reverse b_L and b_N (y and z) vectors in the opposite ("B") strand.
    Matrix<T> basis_reversed; 
    basis_reversed << B.basis.col(0), -B.basis.col(1), -B.basis.col(2);
    const Frame<T> B_reversed = { .basis = basis_reversed, .origin = B.origin };
    return midframe(A, B_reversed);
}


template <typename T>
struct StepParams {
    const T tilt;
    const T roll;
    const T twist;

    const T shift;
    const T slide;
    const T rise;

    const T tiltDeg()  const { return tilt  * 180.0 / PI; }
    const T rollDeg()  const { return roll  * 180.0 / PI; }
    const T twistDeg() const { return twist * 180.0 / PI; }
};

template <typename T>
std::ostream &operator<<(std::ostream &os, StepParams<T> const &f) { 
    return os << "Step params:" << std::endl 
    << "tilt: " << f.tiltDeg()  << "deg" << std::endl 
    << "roll: " << f.rollDeg()  << "deg" << std::endl 
    << "Twist: " << f.twistDeg() << "deg" << std::endl 
    << "shift: " << f.shift << std::endl 
    << "slide: " << f.slide << std::endl 
    << "rise: " << f.rise << std::endl;
}

// Returns step parameters for the two given base-pair frames A and B.
template <typename T>
static inline StepParams<T> stepParams(const Frame<T> &A, const Frame<T> &B) {
    const Matrix<T> Q(rotmatFromFrames(A, B));
    const auto [axis, angle] = axisAngleFromMatrix(Q);

    const Matrix<T> Mb(rotatedHalfway(axis, angle, A.basis));  // mid-basis
    const Vector<T> translation(B.origin - A.origin);

    const Vector<T> x(Mb.col(0));
    const Vector<T> y(Mb.col(1));
    const Vector<T> z(Mb.col(2));

    const StepParams<T> p = {
        .tilt  = angle * axis.dot(x),
        .roll  = angle * axis.dot(y),
        .twist = angle * axis.dot(z),
        .shift = translation.dot(x),
        .slide = translation.dot(y),
        .rise  = translation.dot(z)
    };
    return p;
}

}  // namespace

using DnaParams::Detail::Vector;
using DnaParams::Detail::Matrix;
using DnaParams::Detail::Frame;
using DnaParams::Detail::frameSb;
using DnaParams::Detail::frameBp;
using DnaParams::Detail::StepParams;
using DnaParams::Detail::stepParams;
using DnaParams::Detail::operator<<;

}  // namespace

#endif
