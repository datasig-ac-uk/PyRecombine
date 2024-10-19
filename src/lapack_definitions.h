//
// Created by sam on 19/10/24.
//

#ifndef LAPACK_DEFINITIONS_H
#define LAPACK_DEFINITIONS_H

#include "recombine_common.h"


#ifdef RECOMBINE_FORTAN_LAPACK

extern "C" void DAXPY(integer* N, doublereal* A, doublereal* X, integer* INCX, doublereal* Y, integer* INCY);

extern "C" void DGEMM(char* transa, char* transb, integer* m, integer* n, integer* k, doublereal* alpha, doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal* beta, doublereal* c, integer* ldc);

extern "C" doublereal DASUM(integer* n, doublereal* x, integer* incx);

extern "C" integer IDAMAX(integer* n, doublereal* x, integer* incx);

extern "C" void DGELS(char* transa, integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
                      doublereal* work, integer* lwork, integer* info);

extern "C" void DGELSS(integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
                       doublereal*
                               s,
                       doublereal* rcond, integer* rank, doublereal* work, integer* lwork, integer* info);

extern "C" void DGESVD(char* jobu, char* jobvt, integer* m, integer* n, doublereal* a, integer* lda, doublereal* s, doublereal* u, integer* ldu, doublereal* vt, integer* ldvt, doublereal* work, integer* lwork, integer* info);

extern "C" void DGELSD(integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
                       doublereal* s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork, integer* iwork, integer* info);


#else

// Y=a*X+Y
extern "C" void daxpy_(integer* N, doublereal* A, doublereal* X, integer* INCX, doublereal* Y, integer* INCY);
// returns X.Y
extern "C" doublereal ddot_(integer* N, doublereal* X, integer* INCX, doublereal* Y, integer* INCY);
//minimize 2-norm(| b - A*x |)
extern "C" void dgelsd_(integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
                        doublereal* s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork, integer* iwork, integer* info);
extern "C" void dgels_(char* transa, integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
                       doublereal* work, integer* lwork, integer* info);
extern "C" void dgelss_(integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
                        doublereal* s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork, integer* info);
extern "C" doublereal dnrm2_(integer* n, doublereal* x, integer* incx);
extern "C" void dgesv_(integer* n, integer* nrhs, doublereal* a, integer* lda, integer* ipiv, doublereal* b, integer* ldb, integer* info);
extern "C" void dgesvd_(char* jobu, char* jobvt, integer* m, integer* n, doublereal* a, integer* lda, doublereal* s, doublereal* u, integer* ldu, doublereal* vt, integer* ldvt, doublereal* work, integer* lwork, integer* info);
extern "C" void dcopy_(integer* n, doublereal* x, integer* incx, doublereal* y, integer* incy);
extern "C" void dgemm_(char* transa, char* transb, integer* m, integer* n, integer* k, doublereal* alpha, doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal* beta, doublereal* c, integer* ldc);


#define DGEMM dgemm_
#define DGESVD dgesvd_
#define DGELSD dgelsd_



#endif





#endif //LAPACK_DEFINITIONS_H
