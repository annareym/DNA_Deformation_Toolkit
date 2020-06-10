#include <iostream>
#include <stdexcept>
#include "readpdb.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Expected 1 argument, got " << argc << std::endl;
        exit(1);
    }
    const std::string filename(argv[1]);

    Protein p;
    p.readPDB(filename);

    if (p.chains.size() != 2) {
        throw std::runtime_error("Expected two chains, got " + std::to_string(p.chains.size()));
    }

    auto A = p.chains.at("A");
    auto B = p.chains.at("B");

    if (A.residues.size() != B.residues.size()) {
        throw std::runtime_error("Chain sizes differ: " + std::to_string(A.residues.size()) + \
                                "!=" + std::to_string(B.residues.size()));
    }

    auto itA = A.residues.begin();   // iterator, keys are sorted by default.
    auto itB = B.residues.rbegin();  // reverse iterator.
    for (int i = 0; i < A.residues.size(); i++) {
        auto resA = itA->second;  //  1,  2
        auto resB = itB->second;  // 24, 23

        std::cout << resA.get_C1prime().pdbLine << std::endl;
        std::cout << resA.get_YN1_RN9().pdbLine << std::endl;
        std::cout << resA.get_YC2_RC4().pdbLine << std::endl;

        std::cout << resB.get_C1prime().pdbLine << std::endl;
        std::cout << resB.get_YN1_RN9().pdbLine << std::endl;
        std::cout << resB.get_YC2_RC4().pdbLine << std::endl;

        itA++;
        itB++;
    }
    std::cout << std::endl;
    
    itA = A.residues.begin();
    itB = B.residues.rbegin();
    for (int i = 0; i < A.residues.size(); i++) {
        auto resA = itA->second;  //  1,  2
        auto resB = itB->second;  // 24, 23
        std::cout << resA.get_C1prime().atomNumber << " ";
        std::cout << resA.get_YN1_RN9().atomNumber << " ";
        std::cout << resA.get_YC2_RC4().atomNumber << " ";
        std::cout << resB.get_C1prime().atomNumber << " ";
        std::cout << resB.get_YN1_RN9().atomNumber << " ";
        std::cout << resB.get_YC2_RC4().atomNumber << " ";
        itA++;
        itB++;
    }
    std::cout << std::endl;
}
