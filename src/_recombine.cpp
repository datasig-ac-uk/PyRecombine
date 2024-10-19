// private headers
#include "_recombine.h"
#include "recombine_internal.h"
//#include "TestVec/RdToPowers.h" // CMultiDimensionalBufferHelper
//#include "TestVec/EvaluateAllMonomials.h" //EvaluateAllMonomials::F





void _recombineC(size_t stCubatureDegree, ptrdiff_t dimension, ptrdiff_t no_locations, ptrdiff_t* pno_kept_locations,
                 const void** ppLocationBuffer, double* pdWeightBuffer, size_t* KeptLocations, double* NewWeights)
{
    ptrdiff_t& no_kept_locations = *pno_kept_locations;
    // the required max size of the out buffer needs to be known in advance by the calling routine

    size_t iNoDimensionsToCubature = RdToPowersCubatureDimension(dimension, stCubatureDegree);
    if (0 == no_locations)
    {
        no_kept_locations = iNoDimensionsToCubature;
        return;
    }

    // set up the input structure for conditioning the helper function
    sCMultiDimensionalBufferHelper sConditioning;
    sConditioning.D = stCubatureDegree;
    sConditioning.L = dimension;

    // set up the input structure for data reduction "in"
    sCloud in;

    // chain optional extension information used to condition the data
    in.end = &sConditioning;

    // place the sizes of the buffers and their locations into the structure "in"
    in.NoActiveWeightsLocations = no_locations;
    in.LocationBuf = ppLocationBuffer;
    in.WeightBuf = pdWeightBuffer;

    // set up the output structure for data reduction "out"
    sRCloudInfo out;
    out.end = 0;

    // set the locations of these buffers into the structure "out"
    out.KeptLocations = KeptLocations;
    out.NewWeightBuf = NewWeights;

    // check the sizes of the out buffers
    if (*pno_kept_locations < iNoDimensionsToCubature)
    {
        *pno_kept_locations = 0;
        return;
    }
    // buffers reported big enough
    out.No_KeptLocations = iNoDimensionsToCubature;

    // setup the Recombine Interface data which will join the input and output
    sRecombineInterface data;
    data.end = 0;

    // bind in and out together in data
    data.pInCloud = &in;
    data.pOutCloudInfo = &out;

    // add the degree of the vectors used and the callback function that expands
    // the array of pointers to points into a long buffer of vectors
    data.degree = iNoDimensionsToCubature;

    data.expander = &RdToPowers;

    {
        // CALL THE LIBRARY THAT DOES THE HEAVYLIFTING
        Recombine(&data);
    }
    // recover the information and resize buffers down to the data
    *pno_kept_locations = data.pOutCloudInfo->No_KeptLocations;
}

