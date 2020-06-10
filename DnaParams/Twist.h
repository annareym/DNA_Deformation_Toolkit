#ifndef Twist_h
#define Twist_h

#include <array>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>
#include "DnaParams.h"
#include "Flat.h"

namespace DnaParams {

// Returns total twist for given basepairs.
template <typename T>
static inline T totalTwistDeg(const std::vector<Frame<T> > &bpFrames) {
    T twist = T(0.0);
    for(int i = 0; i < std::size(bpFrames)-1; i++) {
        const StepParams<T> param(stepParams(bpFrames[i], bpFrames[i+1]));
        twist += param.twistDeg();
    }
    return twist;
}

// Copy of totalTwistDeg, with static size.
template <typename T, std::size_t N_BP_FRAMES>
static inline T totalTwistDegSs(const std::array<Frame<T>, N_BP_FRAMES> &bpFrames) {
    T twist = T(0.0);
    for(int i = 0; i < std::size(bpFrames)-1; i++) {
        const StepParams<T> param(stepParams(bpFrames[i], bpFrames[i+1]));
        twist += param.twistDeg();
    }
    return twist;
}

template <typename T>
static inline T totalTwistDegFlat(const Eigen::Matrix<T,Eigen::Dynamic,1> &V) {
    const auto bpFrames = makeBpFrames(V);
    return totalTwistDeg(bpFrames);
}

// Copy of totalTwistDegFlat with static size.
template <typename T, int N_COORDS>
static inline T totalTwistDegFlatSs(const Eigen::Matrix<T, N_COORDS, 1> &V) {
    const auto bpFrames = makeBpFramesSs<T, N_COORDS>(V);
    return totalTwistDegSs(bpFrames);
}

} // namespace
#endif
