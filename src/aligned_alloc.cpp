//
// Created by sam on 19/10/24.
//

#include "aligned_vec.h"


#include <cstdlib>



void *recombine::dtl::aligned_alloc(size_t alignment, size_t size) {
#ifdef WIN32
    return _aligned_malloc(size, alignment);
#else
    return std::aligned_alloc(alignment, size);
#endif
}

void recombine::dtl::aligned_free(void* ptr, size_t size) {
#ifdef WIN32
    _aligned_free(ptr);
#else
    std::free(ptr);
#endif
}