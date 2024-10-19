//
// Created by sam on 18/10/24.
//

#ifndef LINEARALGEBRAREDUCTIONTOOL_H
#define LINEARALGEBRAREDUCTIONTOOL_H



#include "recombine_common.h"
#include "aligned_vec.h"
#include "checked_assign.h"

#include <cassert>

#ifndef REDUCTION_ALGO
#define REDUCTION_ALGO svd
#endif

namespace recombine {

class LinearAlgebraReductionTool {

    VECTORD vdWork;
    VECTORI viWork;

    // counts the number of calls to the linear reduction package
    size_t iNoCallsLinAlg;
    integer iNoCoords;
    integer iNoPoints;
    integer iNoRhs;

public:
    enum MoveMass_type {
        svd,
        simplex
    };


private:
    MoveMass_type MoveMassAlgo;

public:
    LinearAlgebraReductionTool()
        : MoveMassAlgo(REDUCTION_ALGO),
          iNoCoords(1),
          iNoPoints(1),
          iNoRhs(1),
          iNoCallsLinAlg(0)
    {}

    [[nodiscard]]
    integer INoCoords() const
    {
        return iNoCoords;
    }
    const integer& INoCoords(ptrdiff_t val)
    {
        checked_assign(iNoCoords, val);
        return iNoCoords;
    }
    [[nodiscard]]
    integer INoPoints() const
    {
        return iNoPoints;
    }

    const integer& INoPoints(ptrdiff_t val)
    {
        checked_assign(iNoPoints, val);
        return iNoPoints;
    }
    [[nodiscard]]
    size_t INoCallsLinAlg() const
    {
        return iNoCallsLinAlg;
    }

    void MoveMass(VECTORD& eWeights, const VECTORD& ePoints,//aligned_vec<doublereal>& eMassCog,
                  VECTORI& maxset);

private:
    void find_kernel(VECTORD A, integer rowsA, integer lda, VECTORD& K, integer rowsK, integer ldk);

    void MoveMass_svd(VECTORD& eWeights, const VECTORD& ePoints,//aligned_vec<doublereal>& eMassCog,
                      VECTORI& maxset);

    void SharpenWeights(
            VECTORI& minset,
            VECTORI& maxset,
            const VECTORD& ePoints,
            VECTORD& eWeights,
            VECTORD Mcog);

#ifndef NOSIMPLEX
    void MoveMass_simplex(VECTORD& eWeights, const VECTORD& ePoints,//aligned_vec<doublereal>& eMassCog,
                                                       VECTORI& maxset);
#endif
};


}

#endif //LINEARALGEBRAREDUCTIONTOOL_H
