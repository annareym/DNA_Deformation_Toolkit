#ifndef TwistDerivativesAdept_h
#define TwistDerivativesAdept_h

#include <array>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>
#include "DnaParams.h"
#include "FlatTwist.h"

#include <adept.h>

namespace DnaParams {

// Returns twist value and twist derivatives.
static inline std::tuple<double,Eigen::VectorXd> twistDegAndDerivativesAdept(const Eigen::VectorXd &V) {
    adept::Stack stack;
    using adept::adouble;

    Eigen::Matrix<adept::adouble, Eigen::Dynamic, 1> x(V.size());
    for(int i = 0; i < V.size(); i++) {
        x[i] = V[i];
    }
    stack.new_recording();
    adept::adouble y = totalTwistDegFlat(x); 

    y.set_gradient(1.0);
    stack.compute_adjoint();

    Eigen::VectorXd g(V.size());
    for(int i = 0; i < V.size(); i++) {
        g[i] = x[i].get_gradient();
    }
    return {y.value(), g};
}


}  // namespace



#endif
