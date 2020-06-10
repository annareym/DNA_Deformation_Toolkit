// tests-flatten.cpp
#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include "catch.hpp"
#include "DnaParams.h"
#include "readpdb.h"
#include "TwistDerivatives.h"
#include "StretchDerivatives.h"

using namespace DnaParams;

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>
using namespace autodiff;

// Returns squared difference (y-y0)^2
static inline double sqdiff(const double y, const double y0) {
    return std::pow(y-y0, 2);
}

// Returns derivative of squared difference (y-y0)^2
// Given sqdiff(x) = (tw(x) - tw0)^2, then
// sqdiff'(x) = 2 * (tw(x) - tw0) * tw'(x)
template <typename T>
static inline T sqdiffDerivative(const double y, const double y0, const T &y_prime) {
    return 2 * (y-y0) * y_prime;
}

template <typename T>
Eigen::VectorXd fillFlatVector(const Protein &prot, const std::vector<std::pair<int, int> > &basepairIds) {
    Eigen::VectorXd V(std::size(basepairIds) * 18);
    int i = 0;
    for (auto bp: basepairIds) {
        {
            Residue res = prot.chain("A").residue(bp.first);
            Atom C1p = res.get_C1prime();
            V[i++] = C1p.x;
            V[i++] = C1p.y;
            V[i++] = C1p.z;
            Atom N = res.get_YN1_RN9();
            V[i++] = N.x;
            V[i++] = N.y;
            V[i++] = N.z;
            Atom C = res.get_YC2_RC4();
            V[i++] = C.x;
            V[i++] = C.y;
            V[i++] = C.z;
        }
        {
            Residue res = prot.chain("B").residue(bp.second);
            Atom C1p = res.get_C1prime();
            V[i++] = C1p.x; 
            V[i++] = C1p.y;
            V[i++] = C1p.z;
            Atom N = res.get_YN1_RN9();
            V[i++] = N.x;
            V[i++] = N.y;
            V[i++] = N.z;
            Atom C = res.get_YC2_RC4();
            V[i++] = C.x;
            V[i++] = C.y;
            V[i++] = C.z;
        }
    }
    return V;
}

template <typename T, std::size_t N_BP>
Eigen::Matrix<T,N_BP * 18,1> fillFlatVectorSs(const Protein &prot, const std::array<std::pair<int, int>, N_BP > &basepairIds) {
    Eigen::Matrix<T, N_BP * 18, 1> V;
    int i = 0;
    for (auto bp: basepairIds) {
        {
            Residue res = prot.chain("A").residue(bp.first);
            Atom C1p = res.get_C1prime();
            V[i++] = C1p.x;
            V[i++] = C1p.y;
            V[i++] = C1p.z;
            Atom N = res.get_YN1_RN9();
            V[i++] = N.x;
            V[i++] = N.y;
            V[i++] = N.z;
            Atom C = res.get_YC2_RC4();
            V[i++] = C.x;
            V[i++] = C.y;
            V[i++] = C.z;
        }
        {
            Residue res = prot.chain("B").residue(bp.second);
            Atom C1p = res.get_C1prime();
            V[i++] = C1p.x; 
            V[i++] = C1p.y;
            V[i++] = C1p.z;
            Atom N = res.get_YN1_RN9();
            V[i++] = N.x;
            V[i++] = N.y;
            V[i++] = N.z;
            Atom C = res.get_YC2_RC4();
            V[i++] = C.x;
            V[i++] = C.y;
            V[i++] = C.z;
        }
    }
    return V;
}

TEST_CASE( "Flat single-vector total twist is correct") {
    typedef double MyT;
    Protein prot;
    prot.readPDB("1d29.pdb");

    const double tol = 0.001;

    std::array<std::pair<int, int>, 12> basepairIdsArr {{ 
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
    }};
    std::vector<std::pair<int, int> > basepairIds(basepairIdsArr.begin(), basepairIdsArr.end());

    SECTION("Total twist") {
        typedef double MyT;
        Eigen::Matrix<MyT,Eigen::Dynamic,1> V = fillFlatVector<double>(prot, basepairIds);
        MyT twist = totalTwistDegFlat(V);
        REQUIRE(std::abs(twist - 395.71) < tol);
    }

    SECTION("Flat twist one-param") {
        typedef dual MyT;
        const double targetTotalTwistDeg = 405;
        auto f2 = [&](auto x) { return abs(totalTwistDegFlat<MyT>(x)-targetTotalTwistDeg); };

        autodiff::VectorXdual x(fillFlatVector<double>(prot, basepairIds));    

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

    SECTION("Full autodiff") {
        // Use gradient descent to get twist tw(x) to target twist tw0.
        // Difference/discrepancy: d(x) = (tw(x) - tw0)^2.
        // Gradient of difference: d'(x) = 2 * (tw(x) - tw0) * tw'(x)
        const double targetTotalTwistDeg = 405;

        Eigen::VectorXd x(fillFlatVector<double>(prot, basepairIds));

        double t_start = totalTwistDegFlat<double>(x);
        double d = sqdiff(t_start, targetTotalTwistDeg);
        REQUIRE( d > 3 );  // Check that we start "far enough" from 0.

        // Run Gradient Descent.
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        const double alpha = 0.0001;
        for(int i = 0; i < 600; i++) {
            auto [twist, gradTwist] = twistDegAndDerivatives(x);
            d = sqdiff(twist, targetTotalTwistDeg);
            auto d_prime = sqdiffDerivative(twist, targetTotalTwistDeg, gradTwist);
            x -= d_prime * alpha;        
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cerr << "Time difference Dynamic = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

        REQUIRE( d < 1 );  // Should be closer to 0.    
    }

    SECTION("Full autodiff - static size") {
        // Use gradient descent to get twist tw(x) to target twist tw0.
        // Difference/discrepancy: d(x) = (tw(x) - tw0)^2.
        // Gradient of difference: d'(x) = 2 * (tw(x) - tw0) * tw'(x)
        const double targetTotalTwistDeg = 405;

        constexpr std::size_t size = std::size(basepairIdsArr)*18;
        Eigen::Matrix<double, size, 1> x(fillFlatVectorSs<double>(prot, basepairIdsArr));

        double t_start = totalTwistDegFlatSs<double, size>(x);
        double d = sqdiff(t_start, targetTotalTwistDeg);
        REQUIRE( d > 3 );  // Check that we start "far enough" from 0.

        // Run Gradient Descent.
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        const double alpha = 0.0001;
        for(int i = 0; i < 600; i++) {
            auto [twist, gradTwist] = twistDegAndDerivativesSs<size>(x);
            d = sqdiff(twist, targetTotalTwistDeg);
            auto d_prime = sqdiffDerivative(twist, targetTotalTwistDeg, gradTwist);
            x -= d_prime * alpha;
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cerr << "Time difference  Static = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

        REQUIRE( d < 1 );  // Should be closer to 0.    
    }
}


TEST_CASE( "Flat single-vector total stretch is correct") {
    typedef double MyT;
    Protein prot;
    prot.readPDB("1d29.pdb");

    const double tol = 0.001;

    std::array<std::pair<int, int>, 12> basepairIdsArr {{ 
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
    }};
    std::vector<std::pair<int, int> > basepairIds(basepairIdsArr.begin(), basepairIdsArr.end());

    SECTION("Total stretch") {
        typedef double MyT;
        Eigen::Matrix<MyT,Eigen::Dynamic,1> V = fillFlatVector<double>(prot, basepairIds);
        MyT stretch = totalStretchFlat(V);
        //REQUIRE(stretch == tol);
        REQUIRE(std::abs(stretch - 36.6436) < tol);
    }

    SECTION("Full autodiff stretch") {
        // Use gradient descent to get twist tw(x) to target twist tw0.
        // Difference/discrepancy: d(x) = (tw(x) - tw0)^2.
        // Gradient of difference: d'(x) = 2 * (tw(x) - tw0) * tw'(x)
        const double targetStretch = 40;

        Eigen::VectorXd x(fillFlatVector<double>(prot, basepairIds));

        double t_start = totalStretchFlat<double>(x);
        double d = sqdiff(t_start, targetStretch);
        REQUIRE( d > 3 );  // Check that we start "far enough" from 0.

        // Run Gradient Descent.
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        const double alpha = 0.0005;
        for(int i = 0; i < 600; i++) {
            auto [stretch, grad] = stretchDerivatives(x);
            d = sqdiff(stretch, targetStretch);
            auto d_prime = sqdiffDerivative(stretch, targetStretch, grad);
            x -= d_prime * alpha;        
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cerr << "Time difference Dynamic stretch = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

        REQUIRE( d < 1 );  // Should be closer to 0.    
    }

}
