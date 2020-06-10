// tests-twist.cpp
#include "catch.hpp"
#include "DnaParams.h"
#include "readpdb.h"
#include "baseatomcoords.h"
#include "baseatomcoordstopdb.h"

using namespace DnaParams;

TEST_CASE( "Twist is correct") {
    typedef double MyT;
    std::map<int,BaseAtomCoords<MyT> > residueTriplets = makeAtomsMap<MyT>("1d29.pdb");
    const double tol = 0.001;

    SECTION("Simple single base 1") {

        Matrix<MyT> basis; basis <<
            -0.99836 ,    0.044969  ,  -0.035459 ,
             0.056378,      0.88047 ,    -0.47073,
             0.010052,     -0.47196 ,    -0.88156;

        Vector<MyT> origin; origin <<   16.431, 25.659, 26.555;

        BaseAtomCoords<MyT> triplet = residueTriplets.at(1);
        Frame<MyT> sb1 = frameSb(triplet);
        REQUIRE( basis.isApprox(sb1.basis, tol) );
        REQUIRE( origin.isApprox(sb1.origin, tol) );
    }

    SECTION("Simple single base 24") {
        Matrix<MyT> basis; basis <<
            -0.96881  ,    0.016475,     0.24728 ,    
             0.092449 ,   -0.90174 ,     0.42229 ,    
             0.22993  ,    0.43197 ,     0.87208     ;

        Vector<MyT> origin; origin <<     16.186, 25.491, 26.47  ;
        BaseAtomCoords<MyT> triplet = residueTriplets.at(24);
        Frame<MyT> sb = frameSb(triplet);
        REQUIRE( basis.isApprox(sb.basis, tol) );
        REQUIRE( origin.isApprox(sb.origin, tol) );
    }

    SECTION("Simple single base 2") {
        Matrix<MyT> basis; basis <<
            -0.72987  ,    0.68194   , -0.04744     ,
              0.64103 ,     0.65867  ,  -0.39401    ,
             -0.23744 ,    -0.31798  ,  -0.91788    ;

        Vector<MyT> origin; origin <<     16.267, 24.18, 23.105;
        BaseAtomCoords<MyT> triplet = residueTriplets.at(2);
        Frame<MyT> sb = frameSb(triplet);
        REQUIRE( basis.isApprox(sb.basis, tol) );
        REQUIRE( origin.isApprox(sb.origin, tol) );
    }

    SECTION("Simple single base 23") {
        Matrix<MyT> basis; basis <<
            -0.69457  ,   -0.67426   ,   0.25089   ,    
              0.71668 ,    -0.61805  ,    0.32308  ,    
            -0.062785 ,     0.40421  ,    0.91251  ;

        Vector<MyT> origin; origin <<  16.158, 24.413, 23.155 ; 
        BaseAtomCoords<MyT> triplet = residueTriplets.at(23);
        Frame<MyT> sb = frameSb(triplet);
        REQUIRE( basis.isApprox(sb.basis, tol) );
        REQUIRE( origin.isApprox(sb.origin, tol) );
    }

    SECTION("frameBp 1-24") {
        Matrix<MyT> basis; basis <<
              -0.9898  ,   0.013039   ,  -0.14188  ,  
             0.075896  ,    0.89102   ,  -0.44758  ,  
              0.12059  ,   -0.45378   ,  -0.88292  ;

        Vector<MyT> origin; origin <<    16.309, 25.575, 26.513 ; 
        BaseAtomCoords<MyT> triplet1  = residueTriplets.at(1);
        BaseAtomCoords<MyT> triplet24 = residueTriplets.at(24);
        Frame<MyT> sb1  = frameSb(triplet1);
        Frame<MyT> sb24 = frameSb(triplet24);
        Frame<MyT> mf = frameBp(sb1, sb24);

        REQUIRE( (origin - mf.origin).norm() < tol );
        REQUIRE( (basis  - mf.basis ).norm() < tol );
    }
    SECTION("frameBp 2-23") {
        Matrix<MyT> basis; basis <<
         -0.71714,      0.68057,     -0.15011,       
          0.68062,       0.6376,     -0.36086,       
         -0.14988,     -0.36096,     -0.92046;       

        Vector<MyT> origin; origin <<    16.212, 24.297, 23.13 ; 
        BaseAtomCoords<MyT> triplet2  = residueTriplets.at(2);
        BaseAtomCoords<MyT> triplet23 = residueTriplets.at(23);
        Frame<MyT> sb2  = frameSb(triplet2);
        Frame<MyT> sb23 = frameSb(triplet23);
        Frame<MyT> mf = frameBp(sb2, sb23);

        REQUIRE( (origin - mf.origin).norm() < tol );
        REQUIRE( (basis  - mf.basis ).norm() < tol );
    }
    SECTION("Twist step-by-step") {
        Frame<MyT> sb1  = frameSb(residueTriplets.at( 1));
        Frame<MyT> sb24 = frameSb(residueTriplets.at(24));
        Frame<MyT> sb2  = frameSb(residueTriplets.at( 2));
        Frame<MyT> sb23 = frameSb(residueTriplets.at(23));
        Frame<MyT> bp_1_24 = frameBp(sb1, sb24);
        Frame<MyT> bp_2_23 = frameBp(sb2, sb23);

        StepParams<MyT> param(stepParams(bp_1_24, bp_2_23));
        REQUIRE(std::abs(param.twistDeg() - 41.915) < tol);
        REQUIRE(std::abs(param.rise - 3.581) < tol);
    }
    SECTION("Twist all-in-one") {
        const auto bp1 = frameBp(residueTriplets, 1, 24);
        const auto bp2 = frameBp(residueTriplets, 2, 23);
        StepParams<MyT> param(stepParams(bp1, bp2));
        REQUIRE(std::abs(param.twistDeg() - 41.915) < tol);
        REQUIRE(std::abs(param.rise - 3.581) < tol);
    }
    SECTION("Total twist") {
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
        MyT twist = totalTwistDeg<MyT>(residueTriplets, basepairIds);
        REQUIRE(std::abs(twist - 395.71) < tol);
    }

}

TEST_CASE( "Old Twist is correct") {
    typedef double MyT;
    Protein p;// construct a new Protein object
    p.readPDB("1d29.pdb");// read PDB file
    const double tol = 0.001;

    SECTION("Direct reading atoms (old): Simple single base 1") {

        Matrix<MyT> basis; basis <<
            -0.99836 ,    0.044969  ,  -0.035459 ,
             0.056378,      0.88047 ,    -0.47073,
             0.010052,     -0.47196 ,    -0.88156;

        Vector<MyT> origin; origin <<   16.431, 25.659, 26.555;

        SECTION("via atoms") {
            Vector<MyT> C1p = getCoords<MyT>(p, "A", 1, "C1'");
            Vector<MyT>   N = getCoords<MyT>(p, "A", 1, "N1");
            Vector<MyT>   C = getCoords<MyT>(p, "A", 1, "C2");
            Frame<MyT> sb1 = frameSb(C1p, N, C);
            REQUIRE( basis.isApprox(sb1.basis, tol) );
            REQUIRE( origin.isApprox(sb1.origin, tol) );
        }

        SECTION("via Triplets") {
            BaseAtomCoords<MyT> triplet = getTriplet<MyT>(p, "A", 1);
            Frame<MyT> sb1 = frameSb(triplet);
            REQUIRE( basis.isApprox(sb1.basis, tol) );
            REQUIRE( origin.isApprox(sb1.origin, tol) );
        }
    }

    SECTION("Simple single base 24") {
        Matrix<MyT> basis; basis <<
            -0.96881  ,    0.016475,     0.24728 ,    
             0.092449 ,   -0.90174 ,     0.42229 ,    
             0.22993  ,    0.43197 ,     0.87208     ;

        Vector<MyT> origin; origin <<     16.186, 25.491, 26.47  ;
        BaseAtomCoords<MyT> triplet = getTriplet<MyT>(p, "B", 24);
        Frame<MyT> sb = frameSb(triplet);
        REQUIRE( basis.isApprox(sb.basis, tol) );
        REQUIRE( origin.isApprox(sb.origin, tol) );
    }

    SECTION("Simple single base 2") {
        Matrix<MyT> basis; basis <<
            -0.72987  ,    0.68194   , -0.04744     ,
              0.64103 ,     0.65867  ,  -0.39401    ,
             -0.23744 ,    -0.31798  ,  -0.91788    ;

        Vector<MyT> origin; origin <<     16.267, 24.18, 23.105;
        BaseAtomCoords<MyT> triplet = getTriplet<MyT>(p, "A", 2);
        Frame<MyT> sb = frameSb(triplet);
        REQUIRE( basis.isApprox(sb.basis, tol) );
        REQUIRE( origin.isApprox(sb.origin, tol) );
    }

    SECTION("Simple single base 23") {
        Matrix<MyT> basis; basis <<
            -0.69457  ,   -0.67426   ,   0.25089   ,    
              0.71668 ,    -0.61805  ,    0.32308  ,    
            -0.062785 ,     0.40421  ,    0.91251  ;

        Vector<MyT> origin; origin <<  16.158, 24.413, 23.155 ; 
        BaseAtomCoords<MyT> triplet = getTriplet<MyT>(p, "B", 23);
        Frame<MyT> sb = frameSb(triplet);
        REQUIRE( basis.isApprox(sb.basis, tol) );
        REQUIRE( origin.isApprox(sb.origin, tol) );
    }

    SECTION("frameBp 1-24") {
        Matrix<MyT> basis; basis <<
              -0.9898  ,   0.013039   ,  -0.14188  ,  
             0.075896  ,    0.89102   ,  -0.44758  ,  
              0.12059  ,   -0.45378   ,  -0.88292  ;

        Vector<MyT> origin; origin <<    16.309, 25.575, 26.513 ; 
        BaseAtomCoords<MyT> triplet1  = getTriplet<MyT>(p, "A",  1);
        BaseAtomCoords<MyT> triplet24 = getTriplet<MyT>(p, "B", 24);
        Frame<MyT> sb1  = frameSb(triplet1);
        Frame<MyT> sb24 = frameSb(triplet24);
        Frame<MyT> mf = frameBp(sb1, sb24);

        REQUIRE( (origin - mf.origin).norm() < tol );
        REQUIRE( (basis  - mf.basis ).norm() < tol );
    }
    SECTION("frameBp 2-23") {
        Matrix<MyT> basis; basis <<
         -0.71714,      0.68057,     -0.15011,       
          0.68062,       0.6376,     -0.36086,       
         -0.14988,     -0.36096,     -0.92046;       

        Vector<MyT> origin; origin <<    16.212, 24.297, 23.13 ; 
        BaseAtomCoords<MyT> triplet2  = getTriplet<MyT>(p, "A",  2);
        BaseAtomCoords<MyT> triplet23 = getTriplet<MyT>(p, "B", 23);
        Frame<MyT> sb2  = frameSb(triplet2);
        Frame<MyT> sb23 = frameSb(triplet23);
        Frame<MyT> mf = frameBp(sb2, sb23);

        REQUIRE( (origin - mf.origin).norm() < tol );
        REQUIRE( (basis  - mf.basis ).norm() < tol );
    }
    SECTION("Twist") {
        Frame<MyT> sb1  = frameSb(getTriplet<MyT>(p, "A",  1));
        Frame<MyT> sb24 = frameSb(getTriplet<MyT>(p, "B", 24));
        Frame<MyT> sb2  = frameSb(getTriplet<MyT>(p, "A",  2));
        Frame<MyT> sb23 = frameSb(getTriplet<MyT>(p, "B", 23));
        Frame<MyT> bp_1_24 = frameBp(sb1, sb24);
        Frame<MyT> bp_2_23 = frameBp(sb2, sb23);

        StepParams<MyT> param(stepParams(bp_1_24, bp_2_23));
        REQUIRE(std::abs(param.twistDeg() - 41.915) < tol);
        REQUIRE(std::abs(param.rise - 3.581) < tol);
    }

}
