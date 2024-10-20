//
// Created by sam on 19/10/24.
//

#include "aligned_vec.h"


#include <cstdlib>



void *recombine::dtl::aligned_alloc(size_t alignment, size_t size) {
#ifdef WIN32
    return _aligned_malloc(size, alignment);
#elif defined(__linux__) || defined(__APPLE__)
    void * ptr;
    if (posix_memalign(&ptr, alignment, size)) {
        ptr = nullptr;
    }
    return ptr;
#else
    return ::malloc(size);
#endif
}

void recombine::dtl::aligned_free(void* ptr, size_t size) {
#ifdef WIN32
    _aligned_free(ptr);
#elif defined(__linux__) || defined(__APPLE__)
    ::free(ptr);
#else
    return ::free(ptr)
#endif
}