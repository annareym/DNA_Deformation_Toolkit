#ifndef BaseAtomCoordsToPdb_H
#define BaseAtomCoordsToPdb_H

#include "DnaParams.h"
#include "baseatomcoords.h"
#include "readpdb.h"

using namespace DnaParams;

// Returns a vector of coordiantes for a given atom.
template <typename T>
inline Vector<T> getCoords(const Protein& p, const string &chain, const int res, const string &atom) {
    const Atom a = p.chain(chain).residue(res).atom(atom);
    Vector<T> V; V << a.x, a.y, a.z;
    return V;
}

template <typename T>
inline Vector<T> fromAtom(const Atom &a) {
    Vector<T> V; V << a.x, a.y, a.z;
    return V;
}

template <typename T>
inline BaseAtomCoords<T> getTriplet(const Protein& p, const string &chain, const int res) {
    const BaseAtomCoords<T> triplet = {
        .C1p = fromAtom<T>(p.chain(chain).residue(res).get_C1prime()),
        .N   = fromAtom<T>(p.chain(chain).residue(res).get_YN1_RN9()),
        .C   = fromAtom<T>(p.chain(chain).residue(res).get_YC2_RC4()),
    };

    return triplet;
}


template <typename T>
std::map<int,BaseAtomCoords<T> > makeAtomsMap(const Protein& p) {
    std::map<int,BaseAtomCoords<T> >  atomsMap;

    for(auto itChain: p.chains) {
        const std::string chain = itChain.first;
        for(auto itRes: p.chain(chain).residues) {
            const int res = itRes.first;

            // Fail-fast if there is a duplicate residue number.
            if(atomsMap.count(res) > 0) {
               throw "error: residue " + std::to_string(res) + " already exists";
            }

            atomsMap[res] = getTriplet<T>(p, chain, res);
        }
    }

    return atomsMap;
}

template <typename T>
std::map<int,BaseAtomCoords<T> > makeAtomsMap(const std::string filename) {
    Protein p;
    p.readPDB(filename);
    return makeAtomsMap<T>(p);
}

#endif
