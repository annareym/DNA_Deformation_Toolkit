#include "catch.hpp"
#include "DnaParams.h"
#include "readpdb.h"
#include "baseatomcoords.h"
#include "baseatomcoordstopdb.h"

using namespace DnaParams;

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>
using namespace autodiff;

TEST_CASE( "Autodiff gradient descent reduces twist") {

    const std::vector<std::pair<int, int> > basepairIds = { 
        { 1, 24}, 
        { 2, 23}, 
        { 3, 22},
        { 4, 21},
        { 5, 20},
        { 6, 19},
        { 7, 18},
        { 8, 17},
        { 9, 16},
        {10, 15},
        {11, 14},
        {12, 13},
    };

    typedef dual MyT;
    std::map<int,BaseAtomCoords<MyT> > residueTriplets = makeAtomsMap<MyT>("1d29.pdb");
    const auto pair = toFlatAtoms(residueTriplets);

    // using namespace std::placeholders;
    // auto f = std::bind(totalTwistDegFlat<MyT>, _1, pair.second, basepairIds);

    const double targetTotalTwistDeg = 405;
    auto f2 = [&](auto x) { return abs(totalTwistDegFlat<MyT>(x, pair.second, basepairIds)-targetTotalTwistDeg); };

    autodiff::VectorXdual x(pair.first);    

    REQUIRE( f2(x) > 3 );  // Check that we start "far enough" from 0.

    // Run Gradient Descent.
    const double alpha = 0.0001;
    for(int i = 0; i < 50; i++) {
        dual u;
        Eigen::VectorXd g = gradient(f2, wrt(x), at(x), u);
        x -= g * alpha;        
    }
    REQUIRE( f2(x) < 1 );  // Should be closer to 0.
}
