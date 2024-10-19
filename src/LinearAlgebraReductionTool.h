//
// Created by sam on 18/10/24.
//

#ifndef LINEARALGEBRAREDUCTIONTOOL_H
#define LINEARALGEBRAREDUCTIONTOOL_H



#include "recombine_common.h"
#include "aligned_vec.h"

#ifndef REDUCTION_ALGO
#define REDUCTION_ALGO svd
#endif

namespace recombine {

class LinearAlgebraReductionTool {
    ptrdiff_t iNoCoords;
    ptrdiff_t iNoPoints;
    ptrdiff_t iNoRhs;

    using VECTORD = aligned_vec<doublereal>;
    using VECTORI = aligned_vec<integer>;

    VECTORD vdWork;
    VECTORI viWork;

    // counts the number of calls to the linear reduction package
    size_t iNoCallsLinAlg;

public:
    enum MoveMass_type {
        svd,
        simplex
    };

    typedef void (LinearAlgebraReductionTool::*MoveMassFn_t)(VECTORD& eWeights, const VECTORD& ePoints,//aligned_vec<doublereal>& eMassCog,
                                                             VECTORI& maxset);

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

    inline ptrdiff_t INoCoords() const
    {
        return iNoCoords;
    }
    inline const ptrdiff_t& INoCoords(ptrdiff_t val)
    {
        return iNoCoords = val;
    }
    inline ptrdiff_t INoPoints() const
    {
        return iNoPoints;
    }
    inline const ptrdiff_t& INoPoints(ptrdiff_t val)
    {
        return iNoPoints = val;
    }
    inline size_t INoCallsLinAlg() const
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
