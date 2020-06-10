#ifndef Readpdb_h
#define Readpdb_h

#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <stdexcept>

#include "trimstring.h"

using std::string;

static inline std::set<string> split(const std::string &s, const char delim=' ') {
    // https://stackoverflow.com/a/27511119
    std::set<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.insert(std::move(item));
    }
    return elems;    
}

static inline bool contains(const std::set<std::string> &names, const std::string &name) {
    return names.find(name) != names.end();
}

struct Atom {
    float x,y,z;
    int atomNumber;
    std::string pdbLine;
};

struct Residue {
    std::map<string, Atom> atoms;
    std::string name;
    inline Atom atom(const std::string &atomName) const { return atoms.at(atomName); };

    // A,G = purine, C,T = pyrymidine. (R=puRine, Y=pYrimidine).
    inline bool isPurine() const {
        return contains(split("A G DG DG3 DG5 GUA DA DA3 DA5 ADE"), name);
    }

    inline bool isPyrimidine() const { 
        return contains(split("C T DT DT3 DT5 THY DC DC3 DC5 CYT"), name);
    }

    inline Atom get_C1prime() const { return atom("C1'"); };

    inline Atom get_YN1_RN9() const {
        const auto atomName = isPurine() ? "N9" : isPyrimidine() ? "N1" : \
                throw std::runtime_error("Unrecognized residue name " + name);
        return atom(atomName);
    };

    inline Atom get_YC2_RC4() const {
        const auto atomName = isPurine() ? "C4" : isPyrimidine() ? "C2" : \
                throw std::runtime_error("Unrecognized residue name " + name);
        return atom(atomName);
    };
};

struct Chain {
    std::map<int, Residue> residues;
    inline Residue residue(const int resNumber) const { return residues.at(resNumber); };
};

struct Protein {
    std::map<string, Chain> chains;
    
    inline Chain chain(const string &chainName) const { return chains.at(chainName); };
    
    void checkResidues(const string chainName, const string residueName, const int residueNumber) {
        const bool residueExists = chains[chainName].residues.count(residueNumber) > 0;
        if (residueExists) {
            const auto storedResidueName = chains[chainName].residues[residueNumber].name;
            if (storedResidueName != residueName) {
                throw std::runtime_error("Residue number " + std::to_string(residueNumber) + \
                        " got new name " + residueName);
            }
        } else {
            // Residue does not exist, create it to record residue name.
            Residue res = { .name = residueName };
            chains[chainName].residues[residueNumber] = res;
        }
    }

    void checkAtomIsDuplicate(const string chainName, const int residueNumber, const string atomName) {
        bool atomExists = chains[chainName].residues[residueNumber].atoms.count(atomName) > 0;
        if (atomExists) {
            throw std::runtime_error("Atom " + atomName + " already exists in chain " +\
                    chainName + ", residue " + std::to_string(residueNumber));
        }
    }

    void readPDB(const string &filename) {
        std::ifstream file(filename);
        std::string str; 
        while (std::getline(file, str))
        {
            if (str.substr(0,6) != "ATOM  ") {
                continue;
            }

            const int    atomNumber    = std::stoi(str.substr( 6,5));
            const string atomName      = trim_copy(str.substr(11,5));
            const string residueName   = trim_copy(str.substr(17,4));
            const string chainName     = trim_copy(str.substr(21,1));
            const int    residueNumber = std::stoi(str.substr(22,5));
            const float x = std::stof(str.substr(30,8));
            const float y = std::stof(str.substr(38,8));
            const float z = std::stof(str.substr(46,8));

            checkResidues(chainName, residueName, residueNumber);
            checkAtomIsDuplicate(chainName, residueNumber, atomName);

            Atom atom = { x,y,z, atomNumber, str };
            chains[chainName].residues[residueNumber].atoms[atomName] = atom;
        }
    }
};


#endif
