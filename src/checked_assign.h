//
// Created by sam on 19/10/24.
//

#ifndef CHECKED_ASSIGN_H
#define CHECKED_ASSIGN_H

#include <cassert>
#include <limits>
#include <type_traits>


namespace recombine {


template <typename I, typename J>
std::enable_if_t<std::is_integral_v<I> && std::is_integral_v<J>>
checked_assign(I& dst, J src)  {
    assert(src >= std::numeric_limits<I>::min() && src <= std::numeric_limits<I>::max());
    dst = static_cast<I>(src);
}



}

#endif //CHECKED_ASSIGN_H
