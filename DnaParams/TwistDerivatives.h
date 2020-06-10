#ifndef TwistDerivatives_h
#define TwistDerivatives_h

#include <array>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>
#include "DnaParams.h"
#include "Twist.h"

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

namespace DnaParams {



// Returns twist value and twist derivatives.
static inline std::tuple<double,Eigen::VectorXd> twistDegAndDerivatives(const Eigen::VectorXd &V) {
    autodiff::VectorXdual x(V);    
    autodiff::dual u;
    Eigen::VectorXd g = gradient(totalTwistDegFlat<autodiff::dual>, wrt(x), at(x), u);    
    return {val(u), g};
}

// Copy of twistDegAndDerivatives with static size.
template <int N_COORDS>
static inline std::tuple<double,Eigen::Matrix<double,N_COORDS,1>> twistDegAndDerivativesSs(const Eigen::Matrix<double,N_COORDS,1> &V) {
    Eigen::Matrix<autodiff::dual, N_COORDS, 1> x(V);    
    autodiff::dual u;
    Eigen::Matrix<double,N_COORDS,1> g = gradient(totalTwistDegFlatSs<autodiff::dual, N_COORDS>, wrt(x), at(x), u);    
    return {val(u), g};
}

}  // namespace

#endif
