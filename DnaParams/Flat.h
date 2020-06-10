#ifndef Flat_h
#define Flat_h

#include <vector>
#include <Eigen/Core>
#include "DnaParams.h"

namespace DnaParams {

// Returns coordinate frame for a DNA base, where coordinates for C1', N and C atoms are 9 consecutive elements
// in a vector `V` starting at position `i`.
template <typename T>
static inline Frame<T> frameSbFlat(const Eigen::Matrix<T, Eigen::Dynamic, 1> &V, const int i) {
    const Vector<T> C1p = (Vector<T>() << V(i + 0), V(i + 1), V(i + 2)).finished();
    const Vector<T> N   = (Vector<T>() << V(i + 3), V(i + 4), V(i + 5)).finished();
    const Vector<T> C   = (Vector<T>() << V(i + 6), V(i + 7), V(i + 8)).finished();
    return frameSb(C1p, N, C);
}

// Copy of frameSbFlat with static size.
template <typename T, int N_COORDS>
static inline Frame<T> frameSbFlatSs(const Eigen::Matrix<T, N_COORDS, 1> &V, const int i) {
    const Vector<T> C1p = (Vector<T>() << V(i + 0), V(i + 1), V(i + 2)).finished();
    const Vector<T> N   = (Vector<T>() << V(i + 3), V(i + 4), V(i + 5)).finished();
    const Vector<T> C   = (Vector<T>() << V(i + 6), V(i + 7), V(i + 8)).finished();
    return frameSb(C1p, N, C);
}

template <typename T>
static inline std::vector<Frame<T> > makeBpFrames(const Eigen::Matrix<T,Eigen::Dynamic,1> &V) {
    const int atomsPerBase = 9;
    constexpr int atomsPerBasepair = atomsPerBase * 2;

    const int numBasepairs = V.size() / atomsPerBasepair;
    const int remainder    = V.size() % atomsPerBasepair;
    assert(remainder == 0);

    std::vector<Frame<T> > bpFrames(numBasepairs);
    int startIndex = 0;
    for(int bp = 0; bp < numBasepairs; bp++) {

        auto sb1 = frameSbFlat(V, startIndex);
        startIndex += atomsPerBase;
        
        auto sb2 = frameSbFlat(V, startIndex);
        startIndex += atomsPerBase;
        
        bpFrames[bp] = frameBp(sb1, sb2);
    }

    return bpFrames;
}

constexpr int atomsPerBase = 9;
constexpr int atomsPerBasepair = atomsPerBase * 2;

// Copy of makeBpFrames with static size.
template <typename T, int N_COORDS>
static inline std::array<Frame<T>,N_COORDS/atomsPerBasepair> makeBpFramesSs(const Eigen::Matrix<T,N_COORDS,1> &V) {

    constexpr int numBasepairs = N_COORDS / atomsPerBasepair;
    constexpr int remainder    = N_COORDS % atomsPerBasepair;
    static_assert(remainder == 0, "Vector size is not a multiple of number of atoms per base pair");

    std::array<Frame<T>, numBasepairs> bpFrames;
    int startIndex = 0;
    for(int bp = 0; bp < numBasepairs; bp++) {

        auto sb1 = frameSbFlatSs(V, startIndex);
        startIndex += atomsPerBase;
        
        auto sb2 = frameSbFlatSs(V, startIndex);
        startIndex += atomsPerBase;
        
        bpFrames[bp] = frameBp(sb1, sb2);
    }

    return bpFrames;
}
}  // namespace
#endif
