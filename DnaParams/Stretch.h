// Stretch.h
#ifndef Stretch_h
#define Stretch_h

#include <array>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>
#include "DnaParams.h"
#include "Flat.h"

namespace DnaParams {

// Returns total stretch for given basepairs.
template <typename T>
static inline T totalStretch(const std::vector<Frame<T> > &bpFrames) {
    T stretch = T(0.0);
    for(int i = 0; i < std::size(bpFrames)-1; i++) {
        const StepParams<T> param(stepParams(bpFrames[i], bpFrames[i+1]));
        stretch += param.rise;
    }
    return stretch;
}


template <typename T>
static inline T totalStretchFlat(const Eigen::Matrix<T,Eigen::Dynamic,1> &V) {
    const auto bpFrames = makeBpFrames(V);
    return totalStretch(bpFrames);
}


} // namespace
#endif
