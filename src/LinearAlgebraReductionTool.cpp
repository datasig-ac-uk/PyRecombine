//
// Created by sam on 18/10/24.
//

#include "LinearAlgebraReductionTool.h"

#include <algorithm>
#include <cassert>

#include "aligned_vec.h"
#include "reweight.h"
#include "SafeInt3.hpp"
#include "lapack_definitions.h"

using namespace recombine;

namespace {
    // From tjlUtilities.h
    template<class T>
    void fill_index_array(T &index_array) {
        integer ii(0);
        for (auto &i: index_array) {
            i = ii++;
        }
    };
}


void LinearAlgebraReductionTool::MoveMass(VECTORD &eWeights, const VECTORD &ePoints, VECTORI &maxset) {
    switch (MoveMassAlgo) {
#ifndef NOSIMPLEX
    case simplex:
        MoveMass_simplex(eWeights, ePoints, maxset);
        break;
#endif
        case svd:
        default:
            MoveMass_svd(eWeights, ePoints, maxset);
            break;
    }
    iNoCallsLinAlg++;
}

void LinearAlgebraReductionTool::MoveMass_svd(VECTORD &eWeights, const VECTORD &ePoints, VECTORI &maxset) {
    assert(ePoints.size() == iNoCoords * eWeights.size());
    assert(iNoPoints == eWeights.size());

    VECTORI minset;
    maxset.clear();

    VECTORD Mcog(iNoCoords, 0.);
        // compute the conserved quantity
        for (auto i = 0; i < iNoPoints; ++i) {
            for (auto j = 0; j < iNoCoords; ++j) {
                Mcog[j] += ePoints[j + i * iNoCoords] * eWeights[i];
            }
        }

    {
        VECTORD kernel;

        find_kernel(ePoints, iNoCoords, iNoCoords, kernel, iNoPoints, iNoPoints);

        VECTORD weights = eWeights;
        auto P = begin(weights);
        auto K = begin(kernel);
        const integer rP(iNoPoints);
        const integer ldk(iNoPoints); // the stride in K
        const integer cK = SafeInt<integer>(kernel.size() / ldk);
        //SHOW(iNoPoints);
        //SHOW(cK);
        //SHOW(Mcog);
        VECTORI rIDX(iNoPoints), cIDX(cK);
        fill_index_array(rIDX);
        fill_index_array(cIDX);

        auto rIDXb(begin(rIDX)), cIDXb(begin(cIDX));
        const doublereal probability_zero_tolerance(0.);
        reweight(rIDXb, cIDXb, P, K, rP, ldk, cK, probability_zero_tolerance);

        //integer reduced_dimension = rP-cK;
        for (integer i = 0; i < cK; ++i) {
            integer wi = rIDX[i];
            maxset.push_back(wi);
            assert(0. == P[i]);
            eWeights[wi] = 0.;
        }

        for (integer i = cK; i < rP; ++i) {
            // the case where the P is already singular even if the  points are independent
            integer wi = rIDX[i];
            eWeights[wi] = P[i];
            if (0. == eWeights[wi]) {
                maxset.push_back(wi);
            } else {
                minset.push_back(wi);
            }
        }
    }
    //{
    //	// compute the conserved quantity
    //	for (auto i = 0; i < iNoPoints; ++i)
    //	for (auto j = 0; j < iNoCoords; ++j)
    //		Mcog[j] -= ePoints[j + i * iNoCoords] * eWeights[i];
    //	//SHOW(Mcog);
    //}
    // comment out the next line for errors with small numbers
    SharpenWeights(minset, maxset, ePoints, eWeights, std::move(Mcog)); // least squares
}

void LinearAlgebraReductionTool::SharpenWeights(VECTORI &minset, VECTORI &maxset, const VECTORD &ePoints,
                                                VECTORD &eWeights, VECTORD Mcog) {
    for (VECTORI temp_minset; temp_minset.size() < minset.size(); minset.swap(temp_minset)) {
        // INoCoords() >= minset.size() (MxN) since columns are independent
        VECTORD A(iNoCoords * minset.size()), W(minset.size(), 0.), B(Mcog);
        assert(B.size() >= minset.size()); // B has room for the answer
        for (index_integer i = 0; i < minset.size(); ++i) {
            for (integer j = 0; j < iNoCoords; ++j) {
                //W[i] = eWeights[minset[i]],
                // B[j] += ePoints[j + iNoCoords * minset[i]] * eWeights[minset[i]],
                A[j + iNoCoords * i] = ePoints[j + iNoCoords * minset[i]];
            }
        }
        // K was numerical kernel so A has full rank allowing LS

        integer &M(iNoCoords),
                N(SafeInt<integer>(minset.size())),
                NRHS(1),
                &LDA(M),
                &LDB(M), // M >= N
                LWORK,
                INFO(0);

        vdWork.resize(1);
        vdWork[0] = 1.;

#ifndef RECOMBINE_NO_DGELSD
        viWork.resize(1);
        viWork[0] = 1;
        integer RANK;
        LWORK = -1;
        doublereal RCOND(0.);
        std::vector<doublereal> S(N);
        //The singular values of A in decreasing order. The condition number of A in the 2-norm is k2(A) = S(1)/ S(min(m, n))
        if (LWORK == -1) {
            DGELSD(&M,
                   &N,
                   &NRHS,
                   A.data(),
                   &LDA,
                   B.data(),
                   &LDB,
                   S.data(),
                   &RCOND,
                   &RANK,
                   vdWork.data(),
                   &LWORK,
                   viWork.data(),
                   &INFO);
            vdWork.resize(LWORK = (integer) vdWork[0]);
            viWork.resize(viWork[0]);
        }
        DGELSD(&M,
               &N,
               &NRHS,
               A.data(),
               &LDA,
               B.data(),
               &LDB,
               S.data(),
               &RCOND,
               &RANK,
               vdWork.data(),
               &LWORK,
               viWork.data(),
               &INFO);
#else
        LWORK = -1;
        if (LWORK == -1) {
            DGELS("N", &M, &N, &NRHS, &A[0], &LDA, &B[0], &LDB, &WORK[0], &LWORK, &INFO);
            WORK.resize(LWORK = (integer)WORK[0]);
        }
        DGELS("N", &M, &N, &NRHS, &A[0], &LDA, &B[0], &LDB, &WORK[0], &LWORK, &INFO);
#endif

        // B == W approximately
        temp_minset.clear();
        for (index_integer i = 0; i < minset.size(); ++i) {
            if (B[i] <= 0.) {
                eWeights[minset[i]] = 0.;
                maxset.push_back(minset[i]);
            } else {
                eWeights[minset[i]] = B[i];
                temp_minset.push_back(minset[i]);
            }
        }
    }
}

void LinearAlgebraReductionTool::find_kernel(VECTORD A, integer rowsA, integer lda, VECTORD &K, integer rowsK,
                                             integer ldk) {
    integer columnsA((integer) (end(A) - begin(A)) / lda);
    integer ldu(1); // don't return U
    VECTORD u(1);
    VECTORD s(std::min(rowsA, columnsA));
    integer ldvt(columnsA);
    VECTORD vt(ldvt * columnsA);
    integer lwork(-1);
    integer info;

    vdWork.resize(1);
    vdWork[0] = 1.;

    if (lwork == -1) {
        DGESVD((char *) "N",
               (char *) "A",
               &rowsA,
               &columnsA,
               A.data(),
               &lda,
               s.data(),
               u.data(),
               &ldu,
               vt.data(),
               &ldvt,
               vdWork.data(),
               &lwork,
               &info);
        vdWork.resize(lwork = (integer) vdWork[0]);
    }
    DGESVD((char *) "N",
           (char *) "A",
           &rowsA,
           &columnsA,
           A.data(),
           &lda,
           s.data(),
           u.data(),
           &ldu,
           vt.data(),
           &ldvt,
           vdWork.data(),
           &lwork,
           &info);

    auto noNonzeroEV = std::upper_bound(begin(s), end(s), 10e-12, std::greater<>()) - begin(s);
    ////SHOW(s);
    ////SHOW(noNonzeroEV);
    K.resize(ldk * (columnsA - noNonzeroEV));

    for (integer i = noNonzeroEV; i < columnsA; ++i) {
        for (integer j = 0; j < columnsA; ++j) {
            K[j + (i - noNonzeroEV) * ldk] = vt[i + j * ldvt];
        }
    }
}
