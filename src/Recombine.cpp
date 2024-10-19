//
// Created by sam on 18/10/24.
//



#include "recombine_internal.h"


#include <cassert>
#include <map>
#include <memory>
#include <utility>
#include <valarray>
#include <vector>

#include "aligned_vec.h"
#include "TreeBufferHelper.h"
#include "LinearAlgebraReductionTool.h"

using namespace recombine;

namespace internal {

static void ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
                                                    aligned_vec<doublereal>& vdWeightsBuffer,
                                                    std::vector<std::valarray<doublereal>>& vdPointsBuffer);

static void RepackPointBuffer(std::map<ptrdiff_t, ptrdiff_t>& currentroots,
                       std::map<size_t, size_t>& miTreePosition,
                       aligned_vec<doublereal>& weights,
                       aligned_vec<doublereal>& points,
                       size_t pointdimension);

static size_t IdentifyLocationsRemainingAndTheirNewWeights(
        size_t Degree,
        CTreeBufferHelper& bhBufInfo,
        std::map<size_t, size_t>& miTreePosition,
        aligned_vec<doublereal>& vdWeightsBuffer,
        std::vector<std::valarray<doublereal>>& vdPointsBuffer,
        aligned_vec<doublereal>& weights,
        size_t& ICountCalls);

static size_t InsertLeafData(sRecombineInterface& data, std::valarray<doublereal>& vdArrayPointsBuffer,
                             aligned_vec<doublereal>& vdWeightsBuffer);

}


void Recombine(RecombineInterface pInterface) {
    sRecombineInterface& data = *pInterface;

    // expand and insert incoming leaf data into buffers
    std::valarray<doublereal> vdFlatPointsBuffer;
    // InsertLeafData assigns memory: 2 * NPointsIn * data.degree
    // make this a memory mapped file
    aligned_vec<doublereal> vdWeightsBuffer;
    //
    size_t NPointsIn = internal::InsertLeafData(data, vdFlatPointsBuffer, vdWeightsBuffer);
    //
    //
    // Fix the width of DATA (including the leading 1)
    size_t Degree = vdFlatPointsBuffer.size() / vdWeightsBuffer.size();
    assert(data.degree == Degree);

    // reference the locations used for the outgoing data
    size_t& NLocationsKept = (data.pOutCloudInfo)->No_KeptLocations;// number actually returned
    doublereal*& WeightBufOut = (data.pOutCloudInfo)
                                        ->NewWeightBuf;// an external buffer containing the weights of the kept Locations // capacity must be at least iNoDimensionsToCubature + 1
    size_t*& LocationsKept = (data.pOutCloudInfo)
                                            ->KeptLocations;// an external buffer containing the offsets of the kept Locations // capacity must be at least iNoDimensionsToCubature + 1

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // INIT FINISHED //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////// Degree is the max number of of non-degenerate points

    //size_t MaxPoints = Degree + 3;
    //if ( 7 >= NPointsIn )

    // TODO DONT CRASH when rhs is 2*degree
    size_t MaxPoints = 2 * Degree;

    if (1 >= NPointsIn) {
        doublereal* pwOut = WeightBufOut;
        size_t* pl = LocationsKept;
        NLocationsKept = NPointsIn;
        for (size_t iIndex = 0; iIndex < NPointsIn; iIndex++) {
            *(pwOut++) = vdWeightsBuffer[iIndex];
            *(pl++) = iIndex;
        }
    }
    else {

        size_t InitialNoTreesInForest(std::min(MaxPoints, NPointsIn));
        CTreeBufferHelper bhBufInfo(InitialNoTreesInForest, NPointsIn);

        // BAD UNNECCESARY COPY AND MEMORY USE HERE but tricky to remove
        // map buffer to array of val arrays for compatibility reasons
        std::vector<std::valarray<doublereal>>
                vdPointsBuffer(
                        bhBufInfo.end(), std::valarray<doublereal>(std::nan("value not yet assigned"), Degree));
        // populate the leaves

        for (size_t i = 0; i < bhBufInfo.iInitialNoLeaves; i++) {
            std::slice vaSlice(Degree * i, Degree, 1);
            vdPointsBuffer[i] = vdFlatPointsBuffer[vaSlice];
        }

        // now fill the forest using the leafs, leaving exactly iNoTrees roots
        internal::ForestOfWeightedVectorsFromWeightedLeafVectors(bhBufInfo,
                                                                    vdWeightsBuffer,
                                                                    vdPointsBuffer);
        size_t ICountCalls;
        std::map<size_t, size_t> miTreePosition;
        aligned_vec<doublereal> weights;
        //SHOW(NPointsIn);
        ICountCalls = internal::IdentifyLocationsRemainingAndTheirNewWeights(
                Degree,
                bhBufInfo,
                miTreePosition,
                vdWeightsBuffer,
                vdPointsBuffer,
                weights,
                ICountCalls);
        doublereal* pw = WeightBufOut;
        size_t* pl = LocationsKept;
        NLocationsKept = miTreePosition.size();//not weights.size();

        for (auto it = miTreePosition.begin(); it != miTreePosition.end(); ++it) {
            assert(bhBufInfo.isleaf(it->first));
            *(pw++) = weights[it->second];
            *(pl++) = it->first;
        }
        //(size_t iIndex = 0; bhBufInfo.isleaf(iIndex); iIndex++)
    }


}



void internal::ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
                                                    aligned_vec<doublereal>& vdWeightsBuffer,
                                                    std::vector<std::valarray<doublereal>>& vdPointsBuffer)
{
    // create correct initial length and allocate memory for recipient valarrays
    // since slice cannot deduce it
    //// TODO
    //// optimise the OMP so it is faster than the non omp!!
#if 1
    {
        for (size_t iIndex = bhBufInfo.iInitialNoLeaves; iIndex < bhBufInfo.end(); iIndex++) {
            size_t uiLeftParent = bhBufInfo.left(iIndex);
            size_t uiRightParent = bhBufInfo.right(iIndex);
            doublereal left = vdWeightsBuffer[uiLeftParent];
            doublereal right = vdWeightsBuffer[uiRightParent];
            doublereal sum = left + right;
            vdWeightsBuffer[iIndex] = sum;
            std::valarray<doublereal>& dPointsBuffer = vdPointsBuffer[iIndex];
            if (left <= right)
                dPointsBuffer = vdPointsBuffer[uiLeftParent] * (left / sum) + vdPointsBuffer[uiRightParent] * (1 - (left / sum));
            else
                dPointsBuffer = vdPointsBuffer[uiLeftParent] * (1 - (right / sum)) + vdPointsBuffer[uiRightParent] * (right / sum);
        }
    }
#else
    {
        const size_t sz = vdPointsBuffer[0].size(), blocksz(64);
        //#pragma omp parallel for
        for (size_t i = 0; i < sz; i += blocksz)
            for (size_t iIndex = bhBufInfo.iInitialNoLeaves; iIndex < bhBufInfo.end(); iIndex++) {
                std::slice identity(i, std::min(sz, i + blocksz) - i, 1);
                size_t uiLeftParent = bhBufInfo.left(iIndex);
                size_t uiRightParent = bhBufInfo.right(iIndex);
                doublereal left = vdWeightsBuffer[uiLeftParent];
                doublereal right = vdWeightsBuffer[uiRightParent];
                doublereal sum = left + right;
                vdWeightsBuffer[iIndex] = sum;
                std::valarray<doublereal>& dPointsBuffer = vdPointsBuffer[iIndex];
                if (left <= right)
                    dPointsBuffer[identity] = std::valarray<doublereal>(vdPointsBuffer[uiLeftParent][identity]) * (left / sum) + std::valarray<doublereal>(vdPointsBuffer[uiRightParent][identity]) * (1 - (left / sum));
                else
                    dPointsBuffer[identity] = std::valarray<doublereal>(vdPointsBuffer[uiLeftParent][identity]) * (1 - (right / sum)) + std::valarray<doublereal>(vdPointsBuffer[uiRightParent][identity]) * (right / sum);
            }
    }

#endif
}

void internal::RepackPointBuffer(std::map<ptrdiff_t, ptrdiff_t>& currentroots,
                       std::map<size_t, size_t>& miTreePosition, aligned_vec<doublereal>& weights, aligned_vec<doublereal>& points,
                       size_t pointdimension)
{
    std::map<ptrdiff_t, ptrdiff_t> currentrootsnew;
    std::map<size_t, size_t> miTreePositionNew;
    aligned_vec<doublereal> weightsnew(currentroots.size());
    aligned_vec<doublereal> pointsnew(currentroots.size() * pointdimension);

    ptrdiff_t i = 0;
    auto itcurrrts = currentroots.begin();
    for (; itcurrrts != currentroots.end(); ++i, ++itcurrrts) {
        miTreePositionNew[itcurrrts->second] = i;
        currentrootsnew[i] = itcurrrts->second;
        weightsnew[i] = weights[itcurrrts->first];
        for (size_t iM = 0; iM < pointdimension; iM++)
            pointsnew[i * pointdimension + iM] = points[itcurrrts->first * pointdimension + iM];
    }
    points.swap(pointsnew);
    weights.swap(weightsnew);
    currentroots.swap(currentrootsnew);
    miTreePosition.swap(miTreePositionNew);
}

size_t internal::IdentifyLocationsRemainingAndTheirNewWeights(
        size_t Degree,
        CTreeBufferHelper& bhBufInfo,
        std::map<size_t, size_t>& miTreePosition,
        aligned_vec<doublereal>& vdWeightsBuffer,
        std::vector<std::valarray<doublereal>>& vdPointsBuffer,
        aligned_vec<doublereal>& weights,
        size_t& ICountCalls)
{
    /////////////////////////////////////////////////
    //SHOW(vdWeightsBuffer.size());
    //SHOW(vdPointsBuffer.size());

    weights.clear();
    weights.resize(bhBufInfo.iNoTrees);
    // create local buffers
    aligned_vec<doublereal> points(bhBufInfo.iNoTrees * Degree);
    std::map<ptrdiff_t, ptrdiff_t> currentroots;// (bhBufInfo.iNoTrees);
    aligned_vec<ptrdiff_t> maxset;

    bool SomeLinearAlgebraToDo = true;// (bhBufInfo.end() >= bhBufInfo.iNoTrees);
    //assert(SomeLinearAlgebraToDo);

    for (size_t iTreeIndexInFixedBuffer = 0;
         iTreeIndexInFixedBuffer < bhBufInfo.iNoTrees;
         iTreeIndexInFixedBuffer++) {
        ptrdiff_t currentroot = currentroots[iTreeIndexInFixedBuffer] = iTreeIndexInFixedBuffer + bhBufInfo.end() - bhBufInfo.iNoTrees;
        miTreePosition[(size_t)currentroot] = iTreeIndexInFixedBuffer;
        weights[iTreeIndexInFixedBuffer] = vdWeightsBuffer[currentroot];
        for (size_t iM = 0; iM < Degree; iM++)
            points[iTreeIndexInFixedBuffer * Degree + iM] = (vdPointsBuffer[currentroot])[iM];
    }

    //SHOW(miTreePosition.size());
    //SHOW(weights.size());

    recombine::LinearAlgebraReductionTool moLinearAlgebraReductionTool;
    moLinearAlgebraReductionTool.INoCoords(Degree);
    //////////////////// HERE /////////////////////////////////////////
    while (SomeLinearAlgebraToDo) {

        moLinearAlgebraReductionTool.INoPoints(weights.size());
        //moLinearAlgebraReductionTool.INoPoints((ptrdiff_t)bhBufInfo.iNoTrees);
        moLinearAlgebraReductionTool.MoveMass(weights, points, maxset);

        if (maxset.empty()) SomeLinearAlgebraToDo = false;
        while (maxset.size()) {
            size_t togoposition(maxset.back());
            maxset.pop_back();
            miTreePosition.erase(currentroots[togoposition]);
            currentroots.erase(togoposition);
            // if there is at least one non-trivial tree split the last
            // (and so deepest) one to fill vacant slot
            size_t tosplit(miTreePosition.rbegin()->first);
            if (!bhBufInfo.isleaf(tosplit)) {
                size_t tosplitposition = miTreePosition[tosplit];
                miTreePosition.erase(tosplit);
                currentroots.erase(tosplitposition);

                currentroots[togoposition] = bhBufInfo.left(tosplit);
                miTreePosition[bhBufInfo.left(tosplit)] = togoposition;
                weights[togoposition] =
                        weights[tosplitposition] * vdWeightsBuffer[bhBufInfo.left(tosplit)] / vdWeightsBuffer[tosplit];

                currentroots[tosplitposition] = bhBufInfo.right(tosplit);
                miTreePosition[bhBufInfo.right(tosplit)] = tosplitposition;
                weights[tosplitposition] *= vdWeightsBuffer[bhBufInfo.right(tosplit)] / vdWeightsBuffer[tosplit];

                for (size_t iM = 0; iM < Degree; iM++) {
                    points[togoposition * Degree + iM] = (vdPointsBuffer[bhBufInfo.left(tosplit)])[iM];
                    points[tosplitposition * Degree + iM] = (vdPointsBuffer[bhBufInfo.right(tosplit)])[iM];
                }
            }
        }

        RepackPointBuffer(currentroots, miTreePosition, weights, points, Degree);
        ICountCalls = moLinearAlgebraReductionTool.INoCallsLinAlg();
        //SHOW(ICountCalls);
    }

    return ICountCalls;
}

size_t internal::InsertLeafData(sRecombineInterface& data, std::valarray<doublereal>& vdArrayPointsBuffer,
                             aligned_vec<doublereal>& vdWeightsBuffer)
{
    void*& LocationBufIn = (data.pInCloud)->LocationBuf;
    size_t NPointsIn = (data.pInCloud)->NoActiveWeightsLocations;
    vdArrayPointsBuffer.resize(2 * NPointsIn * data.degree, std::nan("a non number"));
    vdWeightsBuffer.resize(2 * NPointsIn, std::nan("a non number"));

    // Buffers large enough for any encompassing tree + 1 unused

    auto PointsToVectorDoubles = &(*data.expander);
    sCConditionedBufferHelper arg3;
    arg3.NoPointsToBeProcessed = NPointsIn;
    arg3.SmallestReducibleSetSize = data.degree + 1;//legacy reasons
    arg3.pvCConditioning = (*data.pInCloud).end;
    PointsToVectorDoubles(LocationBufIn, &vdArrayPointsBuffer[0], &arg3);

    doublereal* WeightBufIn = (data.pInCloud)->WeightBuf;
    std::copy(WeightBufIn, WeightBufIn + NPointsIn, vdWeightsBuffer.begin());

    return NPointsIn;
}
