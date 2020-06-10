#include "catch.hpp"
#include "readpdb.h"

using Catch::Matchers::Approx;

TEST_CASE( "PDB reader can read atoms") {

    Protein p;
    p.readPDB("1d29.pdb");

    std::vector<float> xyzExpected = {20.037,  33.839,  25.487};
    Atom a = p.chain("A").residue(1).atom("O5'");
    std::vector<float> xyzActual = {a.x, a.y, a.z};

    REQUIRE_THAT(xyzExpected, Approx(xyzActual));
}
