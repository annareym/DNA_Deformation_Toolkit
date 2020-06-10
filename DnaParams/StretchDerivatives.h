#ifndef StretchDerivatives_h
#define StretchDerivatives_h

#include <array>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>
#include "DnaParams.h"
#include "Stretch.h"

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

namespace DnaParams {

// Returns Stretch value and derivatives.
static inline std::tuple<double,Eigen::VectorXd> stretchDerivatives(const Eigen::VectorXd &V) {
    autodiff::VectorXdual x(V);    
    autodiff::dual u;
    Eigen::VectorXd g = gradient(totalStretchFlat<autodiff::dual>, wrt(x), at(x), u);    
    return {val(u), g};
}

}  // namespace

#endif
