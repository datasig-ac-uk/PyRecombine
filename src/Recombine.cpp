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

using TreePositionMap = std::map<size_t, size_t>;
using CurrentRootsMap = std::map<integer, ptrdiff_t>;

using PointType = std::valarray<doublereal>;
using PointsBuffer = std::vector<PointType>;


namespace internal {


void ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
                                                    VECTORD& vdWeightsBuffer,
                                                    PointsBuffer& vdPointsBuffer);

void RepackPointBuffer(CurrentRootsMap &currentroots, TreePositionMap &miTreePosition,
                       VECTORD &weights, VECTORD &points, size_t pointdimension);

size_t IdentifyLocationsRemainingAndTheirNewWeights(size_t Degree, CTreeBufferHelper &bhBufInfo,
                                                    TreePositionMap &miTreePosition,
                                                    VECTORD &vdWeightsBuffer,
                                                    PointsBuffer &vdPointsBuffer,
                                                    VECTORD &weights, size_t &ICountCalls);

size_t InsertLeafData(sRecombineInterface& data, PointType& vdArrayPointsBuffer,
                             VECTORD& vdWeightsBuffer);



}






void Recombine(RecombineInterface pInterface) {
  // unpack the void pointer
    sRecombineInterface& data = *pInterface;

    // expand and insert incoming leaf data into buffers
    PointType vdFlatPointsBuffer;
    // InsertLeafData assigns memory: 2 * NPointsIn * data.degree
    // make this a memory mapped file
    VECTORD vdWeightsBuffer;
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
       PointsBuffer
                vdPointsBuffer(
                        bhBufInfo.end(), PointType(std::nan("value not yet assigned"), Degree));
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
        TreePositionMap miTreePosition;
        VECTORD weights;
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

        for (auto & it : miTreePosition) {
            assert(bhBufInfo.isleaf(it->first));
            *(pw++) = weights[it.second];
            *(pl++) = it.first;
        }
        //(size_t iIndex = 0; bhBufInfo.isleaf(iIndex); iIndex++)
    }


}



void internal::ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
                                                    VECTORD& vdWeightsBuffer,
                                                    PointsBuffer& vdPointsBuffer)
{
    // create correct initial length and allocate memory for recipient valarrays
    // since slice cannot deduce it
    //// TODO
    //// optimise the OMP so it is faster than the non omp!!
#if 1
    {
        for (auto iIndex = bhBufInfo.iInitialNoLeaves; iIndex < bhBufInfo.end(); iIndex++) {
            auto uiLeftParent = bhBufInfo.left(iIndex);
            auto uiRightParent = bhBufInfo.right(iIndex);

            auto left = vdWeightsBuffer[uiLeftParent];
            auto right = vdWeightsBuffer[uiRightParent];

            auto sum = left + right;
            vdWeightsBuffer[iIndex] = sum;

            auto& dPointsBuffer = vdPointsBuffer[iIndex];
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
            for (index_integer iIndex = bhBufInfo.iInitialNoLeaves; iIndex < bhBufInfo.end(); iIndex++) {
                std::slice identity(i, std::min(sz, i + blocksz) - i, 1);
                index_integer uiLeftParent = bhBufInfo.left(iIndex);
                index_integer uiRightParent = bhBufInfo.right(iIndex);
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

void internal::RepackPointBuffer(CurrentRootsMap& currentroots,
                       TreePositionMap& miTreePosition, VECTORD& weights, VECTORD& points,
                       size_t pointdimension)
{
    CurrentRootsMap currentrootsnew;
    TreePositionMap miTreePositionNew;
    VECTORD weightsnew(currentroots.size());
    VECTORD pointsnew(currentroots.size() * pointdimension);

    integer i = 0;
    auto itcurrrts = currentroots.begin();
    for (; itcurrrts != currentroots.end(); ++i, ++itcurrrts) {
        miTreePositionNew[itcurrrts->second] = i;
        currentrootsnew[i] = itcurrrts->second;
        weightsnew[i] = weights[itcurrrts->first];

        for (index_integer iM = 0; iM < pointdimension; iM++)
            pointsnew[i * pointdimension + iM] = points[itcurrrts->first * pointdimension + iM];
    }
    points.swap(pointsnew);
    weights.swap(weightsnew);
    currentroots.swap(currentrootsnew);
    miTreePosition.swap(miTreePositionNew);
}

size_t internal::IdentifyLocationsRemainingAndTheirNewWeights(size_t Degree, CTreeBufferHelper &bhBufInfo,
                                                              TreePositionMap &miTreePosition,
                                                              VECTORD &vdWeightsBuffer,
                                                              PointsBuffer &vdPointsBuffer,
                                                              VECTORD &weights, size_t &ICountCalls) {
     /////////////////////////////////////////////////
    //SHOW(vdWeightsBuffer.size());
    //SHOW(vdPointsBuffer.size());

    weights.clear();
    weights.resize(bhBufInfo.iNoTrees);
    // create local buffers
    VECTORD points(bhBufInfo.iNoTrees * Degree);
    CurrentRootsMap currentroots;// (bhBufInfo.iNoTrees);
    VECTORI maxset;

    bool SomeLinearAlgebraToDo = true;// (bhBufInfo.end() >= bhBufInfo.iNoTrees);
    //assert(SomeLinearAlgebraToDo);

    for (integer iTreeIndexInFixedBuffer = 0;
         iTreeIndexInFixedBuffer < bhBufInfo.iNoTrees;
         iTreeIndexInFixedBuffer++) {
        auto currentroot = currentroots[iTreeIndexInFixedBuffer] = iTreeIndexInFixedBuffer + bhBufInfo.end() - bhBufInfo.iNoTrees;
        miTreePosition[(index_integer)currentroot] = iTreeIndexInFixedBuffer;
        weights[iTreeIndexInFixedBuffer] = vdWeightsBuffer[currentroot];

        for (index_integer iM = 0; iM < Degree; iM++)
            points[iTreeIndexInFixedBuffer * Degree + iM] = (vdPointsBuffer[currentroot])[iM];
    }

    //SHOW(miTreePosition.size());
    //SHOW(weights.size());

    integer tosplitposition, togoposition;

    recombine::LinearAlgebraReductionTool moLinearAlgebraReductionTool;
    moLinearAlgebraReductionTool.INoCoords(Degree);
    //////////////////// HERE /////////////////////////////////////////
    while (SomeLinearAlgebraToDo) {

        moLinearAlgebraReductionTool.INoPoints(weights.size());
        //moLinearAlgebraReductionTool.INoPoints((integer)bhBufInfo.iNoTrees);
        moLinearAlgebraReductionTool.MoveMass(weights, points, maxset);

        if (maxset.empty()) SomeLinearAlgebraToDo = false;
        while (maxset.size()) {
            checked_assign(togoposition, maxset.back());
            maxset.pop_back();
            miTreePosition.erase(currentroots[togoposition]);
            currentroots.erase(togoposition);
            // if there is at least one non-trivial tree split the last
            // (and so deepest) one to fill vacant slot
            index_integer tosplit(miTreePosition.rbegin()->first);
            if (!bhBufInfo.isleaf(tosplit)) {
                checked_assign(tosplitposition, miTreePosition[tosplit]);
                miTreePosition.erase(tosplit);
                currentroots.erase(tosplitposition);

                currentroots[togoposition] = bhBufInfo.left(tosplit);
                miTreePosition[bhBufInfo.left(tosplit)] = togoposition;
                weights[togoposition] =
                        weights[tosplitposition] * vdWeightsBuffer[bhBufInfo.left(tosplit)] / vdWeightsBuffer[tosplit];

                currentroots[tosplitposition] = bhBufInfo.right(tosplit);
                miTreePosition[bhBufInfo.right(tosplit)] = tosplitposition;
                weights[tosplitposition] *= vdWeightsBuffer[bhBufInfo.right(tosplit)] / vdWeightsBuffer[tosplit];

                for (index_integer iM = 0; iM < Degree; iM++) {
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

size_t internal::InsertLeafData(sRecombineInterface &data, PointType &vdArrayPointsBuffer,
                                       VECTORD &vdWeightsBuffer)
{
    void*& LocationBufIn = (data.pInCloud)->LocationBuf;
    index_integer NPointsIn = (data.pInCloud)->NoActiveWeightsLocations;
    vdArrayPointsBuffer.resize(2 * NPointsIn * data.degree, std::nan("a non number"));
    vdWeightsBuffer.resize(2 * NPointsIn, std::nan("a non number"));

    // Buffers large enough for any encompassing tree + 1 unused

    auto PointsToVectorDoubles = data.expander;
    sCConditionedBufferHelper arg3;
    arg3.NoPointsToBeProcessed = NPointsIn;
    arg3.SmallestReducibleSetSize = data.degree + 1;//legacy reasons
    arg3.pvCConditioning = (*data.pInCloud).end;
    PointsToVectorDoubles(LocationBufIn, &vdArrayPointsBuffer[0], &arg3);

    doublereal* WeightBufIn = (data.pInCloud)->WeightBuf;
    std::copy(WeightBufIn, WeightBufIn + NPointsIn, vdWeightsBuffer.begin());

    return NPointsIn;
}
