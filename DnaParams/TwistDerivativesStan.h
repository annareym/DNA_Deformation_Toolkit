#ifndef TwistDerivativesAdept_h
#define TwistDerivativesAdept_h

#include <array>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>
#include "DnaParams.h"
#include "FlatTwist.h"

#include <stan/math.hpp>

namespace DnaParams {

// Returns twist value and twist derivatives.
static inline std::tuple<double,Eigen::VectorXd> twistDegAndDerivativesStan(const Eigen::VectorXd &V) {


    Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x(V.size());
    for(int i = 0; i < V.size(); i++) {
        x[i] = V[i];
    }
    stan::math::var y = totalTwistDegFlat(x);

    double ans = v.val();

    y.grad();

    Eigen::VectorXd g(V.size());
    y.grad(x, g);

    return {ans, g};
}


}  // namespace



#endif
