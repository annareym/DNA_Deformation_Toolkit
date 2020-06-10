#ifndef BaseAtomCoords_H
#define BaseAtomCoords_H

#include <Eigen/Core>
#include "DnaParams.h"
#include "TwistDerivatives.h"

using namespace DnaParams;

template <typename T>
struct BaseAtomCoords {
    Vector<T> C1p;
    Vector<T> N;
    Vector<T> C;
};

template <typename T>
static std::ostream &operator<<(std::ostream &os, BaseAtomCoords<T> const &f) { 
    return os << "Atoms:" << std::endl 
    << "C1p " << f.C1p << std::endl 
    << "N " << f.N << std::endl 
    << "C " << f.C << std::endl;
}

template <typename T>
static inline Frame<T> frameSb(const BaseAtomCoords<T> &triplet) {
    return frameSb(triplet.C1p, triplet.N, triplet.C);
}

template <typename T>
static inline Frame<T> frameBp(const std::map<int,BaseAtomCoords<T> > &atomTriplets, const int a, const int b) {
    const Frame<T> cf1 = frameSb(atomTriplets.at(a));
    const Frame<T> cf2 = frameSb(atomTriplets.at(b));
    return frameBp(cf1, cf2);
}

template <typename T>
static inline T totalTwistDeg(const std::map<int,BaseAtomCoords<T> > &atomTriplets, const std::vector<std::pair<int, int> > &basepairIds) {
    std::vector<Frame<T> > bpFrames(basepairIds.size());
    for(int i = 0; i < basepairIds.size(); i++) {
        const auto bpi = basepairIds[i];
        bpFrames[i] = frameBp(atomTriplets, bpi.first, bpi.second);
    }

    return totalTwistDeg(bpFrames);
}


template <typename T>
static std::pair<Eigen::Matrix<T,Eigen::Dynamic,1>, std::map<int, int> > toFlatAtoms(const std::map<int,BaseAtomCoords<T> > &atomTriplets) {
    Eigen::Matrix<T,Eigen::Dynamic,1> V(atomTriplets.size() * 9);
    std::map<int, int> reference;
    int i = 0;
    for(auto itAtoms: atomTriplets) {
        reference[itAtoms.first] = i;
        auto t = itAtoms.second;  // triplet.
        V(i++) = t.C1p(0);
        V(i++) = t.C1p(1);
        V(i++) = t.C1p(2);
        V(i++) = t.N(0);
        V(i++) = t.N(1);
        V(i++) = t.N(2);
        V(i++) = t.C(0);
        V(i++) = t.C(1);
        V(i++) = t.C(2);
    }
    return std::pair<Eigen::Matrix<T,Eigen::Dynamic,1>, std::map<int, int> >(V, reference);
}

template <typename T>
static std::map<int,BaseAtomCoords<T> > fromFlatAtoms(const Eigen::Matrix<T,Eigen::Dynamic,1> &V, const std::map<int,int> &reference) {
    std::map<int,BaseAtomCoords<T> > atomTriplets;
    for(auto it: reference) {
        const int i = it.second;  // map contains starting index (out of 9).
        const BaseAtomCoords<T> triplet = {
            .C1p = (Vector<T>() << V(i+0), V(i+1), V(i+2)).finished(),
            .N   = (Vector<T>() << V(i+3), V(i+4), V(i+5)).finished(),
            .C   = (Vector<T>() << V(i+6), V(i+7), V(i+8)).finished(),
        };
        int resNum = it.first;
        atomTriplets[resNum] = triplet;
    }
    return atomTriplets;
}


template <typename T>
static inline T totalTwistDegFlat(const Eigen::Matrix<T,Eigen::Dynamic,1> &V, const std::map<int,int> &reference, const std::vector<std::pair<int, int> > &basepairIds) {
    const auto atomTriplets = fromFlatAtoms(V, reference);
    return totalTwistDeg(atomTriplets, basepairIds);
}


#endif
