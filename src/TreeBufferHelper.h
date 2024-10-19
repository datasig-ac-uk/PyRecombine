//
// Created by sam on 18/10/24.
//

#ifndef TREEBUFFERHELPER_H
#define TREEBUFFERHELPER_H

#include <algorithm>
#include <cassert>
#include <cstdint>

#include "recombine_common.h"

namespace recombine {

struct CTreeBufferHelper {
    // the number of trees in the initial forest
    ptrdiff_t iNoTrees;
    // the number of leaves in the initial forest
    ptrdiff_t iInitialNoLeaves;
    // vdBuffer[iIndex +  iNoPointsToBeprocessed]
    // = vdBuffer[ 2 * iIndex ] + vdBuffer[ 2 * iIndex + 1 ] ;

    CTreeBufferHelper(ptrdiff_t SmallestReducibleSetSize, ptrdiff_t NoPointsToBeprocessed) :
        iNoTrees(SmallestReducibleSetSize), iInitialNoLeaves(NoPointsToBeprocessed) {
        assert(iInitialNoLeaves >= iNoTrees && iNoTrees > 0);
    }

    [[nodiscard]] bool isleaf(ptrdiff_t node) const noexcept { return node < iInitialNoLeaves && node >= 0; }

    [[nodiscard]] ptrdiff_t end() const noexcept { return 2 * iInitialNoLeaves - iNoTrees; }

    [[nodiscard]] bool isnode(ptrdiff_t node) const noexcept { return node >= 0 && node < end(); }

    [[nodiscard]] ptrdiff_t parent(ptrdiff_t node) const noexcept {
        assert(isnode(node));
        return std::min(iInitialNoLeaves + (node / 2), end());
    }

    [[nodiscard]] bool isroot(ptrdiff_t node) const {
        assert(isnode(node));
        return parent(node) == end();
    }

    [[nodiscard]] ptrdiff_t left(ptrdiff_t node) const {
        assert(isnode(node));
        return (node - iInitialNoLeaves) * 2;
    }

    [[nodiscard]] ptrdiff_t right(ptrdiff_t node) const { return left(node) < 0 ? -1 : left(node) + 1; }
};


} // namespace recombine


#endif // TREEBUFFERHELPER_H
