// tests-flatten.cpp
#include <iostream>
#include "catch.hpp"
#include "DnaParams.h"
#include "readpdb.h"
#include "baseatomcoords.h"
#include "baseatomcoordstopdb.h"
#include "TwistDerivatives.h"

using namespace DnaParams;

TEST_CASE( "Flattening is correct") {
    typedef double MyT;
    std::map<int,BaseAtomCoords<MyT> > residueTriplets = makeAtomsMap<MyT>("1d29.pdb");
    const double tol = 0.001;

    std::vector<std::pair<int, int> > basepairIds = { 
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

    SECTION("Simple single base 1") {
        const auto pair = toFlatAtoms(residueTriplets);
        const auto flat1 = pair.first;
        const auto ref1 = pair.second;
        const auto residueTriplets2 = fromFlatAtoms(flat1, ref1);
        const auto pair2 = toFlatAtoms(residueTriplets2);
        const auto flat2 = pair2.first;
        REQUIRE( flat1 == flat2);
    }

    SECTION("Total twist") {
        typedef double MyT;
        std::map<int,BaseAtomCoords<MyT> > residueTriplets = makeAtomsMap<MyT>("1d29.pdb");
        MyT twist = totalTwistDeg(residueTriplets, basepairIds);
        REQUIRE(std::abs(twist - 395.71) < tol);
    }


    SECTION("Flat twist two-params") {
        typedef double MyT;
        std::map<int,BaseAtomCoords<MyT> > residueTriplets = makeAtomsMap<MyT>("1d29.pdb");
        auto pair = toFlatAtoms(residueTriplets);

        Eigen::VectorXd x(pair.first);
        auto reference = pair.second;

        MyT twist = totalTwistDegFlat(x, reference, basepairIds);

        REQUIRE(std::abs(twist - 395.71) < tol);
    }

    SECTION("Flat twist one-param") {
        typedef double MyT;
        std::map<int,BaseAtomCoords<MyT> > residueTriplets = makeAtomsMap<MyT>("1d29.pdb");
        auto pair = toFlatAtoms(residueTriplets);
        Eigen::VectorXd x(pair.first);
        auto reference = pair.second;

        using namespace std::placeholders;

        auto f = [&](auto x) { return totalTwistDegFlat<MyT>(x, reference, basepairIds); };

        MyT twist = f(x);

        REQUIRE(std::abs(twist - 395.71) < tol);
    }
}
