/**
* This is the combination and clean up of recombine.h, BufferConstructor.h
* from the original recombine library without all the garbage that is just
* complicating the repository.
*/


#ifndef RECOMBINE_INTERNAL_H
#define RECOMBINE_INTERNAL_H

#include "recombine_common.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef double doublereal;

typedef void (*expander_fn)(void*, doublereal*, void*);


struct sRecombineInterface {
    struct sCloud* pInCloud;
    struct sRCloudInfo* pOutCloudInfo;
    size_t degree;
    expander_fn expander;
    void* end;
};

typedef struct sRecombineInterface* RecombineInterface;

struct sCloud {
    size_t NoActiveWeightsLocations;
    doublereal* WeightBuf;
    void* LocationBuf;
    void* end;
};

// a C structure for the returned information used for modifying the cloud
struct sRCloudInfo {
    size_t No_KeptLocations;// number actually returned
    doublereal* NewWeightBuf;      // a buffer containing the weights of the kept Locations // capacity must be at least degree + 1
    size_t* KeptLocations;  // a buffer containing the offsets of the kept Locations // capacity must be at least degree + 1
    void* end;
};

struct sCConditionedBufferHelper {
    size_t SmallestReducibleSetSize;//Target Point Dimension + 1
    size_t NoPointsToBeProcessed;
    void* pvCConditioning;
};

typedef struct sCConditionedBufferHelper* CBufferHelper;

struct sCMultiDimensionalBufferHelper {
    size_t L;
    size_t D;
};

enum ProductType {
    Prods2 = 0,
    Prods_test = 1,
    Prods_nonrecursive3 = 2,
    Prods_nonrecursive2 = 3,
    Prods_nonrecursive = 4,
    Prods_wei1 = 5,
    Prods_cheb = 6,
    Prods = 7
};


void Recombine(RecombineInterface pInterface);

// void Compare(RecombineInterface pInterface);


void RdToPowers(void* pIn, doublereal* pOut, void* vpCBufferHelper);
size_t RdToPowersCubatureDimension(size_t stDimension, size_t stCubatureDegree);


#ifdef __cplusplus
}
#endif



#endif //RECOMBINE_INTERNAL_H
