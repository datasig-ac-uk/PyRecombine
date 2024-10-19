//
// Created by sam on 18/10/24.
//

#ifndef RECOMBINE_COMMON_H
#define RECOMBINE_COMMON_H

#include <stddef.h>
#include <stdint.h>

#ifndef RECOMBINE_INT_SIZE
#define RECOMBINE_INT_SIZE 8
#endif

#if INTPTR_MAX == INT64_MAX
#if RECOMBINE_INT_SIZE == 8
typedef ptrdiff_t integer;
typedef size_t index_integer;
#elif RECOMBINE_INT_SIZE == 4
typedef int32_t integer;
typedef uint32_t index_integer;
#else
#error "Valid choices for integer are 4 bytes or 8 bytes"
#endif
#else
typedef ptrdiff_t integer;
typedef size_t index_integer;
#endif

typedef double doublereal;




#endif //RECOMBINE_COMMON_H
