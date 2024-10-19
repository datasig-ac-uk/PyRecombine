//
// Created by sam on 18/10/24.
//

#ifndef ALIGNED_VEC_H
#define ALIGNED_VEC_H

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>

#include "recombine_common.h"

namespace recombine {


namespace dtl {

template <typename I>
constexpr bool is_power_2(I val) noexcept {
    return (val & (val - 1)) == 0;
}


// This is partly based on aligned_alloc.h from the original repository,
// but more from https://stackoverflow.com/a/24361178

template <typename Ty, size_t Alignment>
class AlignedAllocator
{
    static_assert(is_power_2(Alignment), "Aligned allocator must be a power of 2");
public:
    static constexpr size_t alignment = std::max(Alignment, alignof(Ty));

    typedef Ty value_type;
    typedef Ty* pointer;
    typedef const Ty* const_pointer;
    typedef Ty& reference;
    typedef const Ty& const_reference;

    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    typedef std::true_type propagate_on_container_move_assignment;

    template <typename Other>
    struct rebind {
        typedef AlignedAllocator<Other, Alignment> other;
    };

    constexpr AlignedAllocator() noexcept {}

    template <typename Other>
    constexpr AlignedAllocator(const AlignedAllocator<Other, alignment>&) noexcept = delete;

    [[nodiscard]] constexpr size_type max_size() const noexcept {
        return (~static_cast<size_type>(0) - alignment) / sizeof(Ty);
    }

    pointer address(reference val) const noexcept {
        return std::addressof(val);
    }

    const_pointer address(const_reference val) const noexcept {
        return std::addressof(val);
    }

    pointer allocate(size_type n, const void* hint=nullptr) {
        void* ptr =  std::aligned_alloc(alignment, n * sizeof(Ty));
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }
        return static_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type N) {
        (void) N;
        std::free(p);
    }

    template <typename Oy, typename... Args>
    void construct(Oy* p, Args&&... args) {
        ::new(reinterpret_cast<void*>(p)) Oy(std::forward<Args>(args)...);
    }

    void destroy(pointer p) {
        p->~Ty();
    }
};

template <typename Ty, size_t Alignment>
class AlignedAllocator<const Ty, Alignment>
{
    static_assert(is_power_2(Alignment), "Aligned allocator must be a power of 2");
public:
    static constexpr size_t alignment = std::max(Alignment, alignof(Ty));

    typedef Ty const value_type;
    typedef const Ty* pointer;
    typedef const Ty* const_pointer;
    typedef const Ty& reference;
    typedef const Ty& const_reference;

    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    typedef std::true_type propagate_on_container_move_assignment;

    template <typename Other>
    struct rebind {
        typedef AlignedAllocator<Other, Alignment> other;
    };

    constexpr AlignedAllocator() noexcept {}

    template <typename Other>
    constexpr AlignedAllocator(const AlignedAllocator<Other, alignment>&) noexcept = delete;

    [[nodiscard]] constexpr size_type max_size() const noexcept {
        return (~static_cast<size_type>(0) - alignment) / sizeof(Ty);
    }

    pointer address(reference val) const noexcept {
        return std::addressof(val);
    }

    pointer allocate(size_type n, const void* hint=nullptr) {
        void* ptr = std::aligned_alloc(alignment, n * sizeof(Ty));
        if (ptr == nullptr) {
            throw std::bad_alloc {};
        }
        return static_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type N) {
        (void) N;
        std::free(p);
    }

    template <typename Oy, typename... Args>
    void construct(Oy* p, Args&&... args) {
        ::new(reinterpret_cast<void*>(p)) Oy(std::forward<Args>(args)...);
    }

    void destroy(pointer p) {
        p->~Ty();
    }
};

template <size_t Alignment>
class AlignedAllocator<void, Alignment> {
    static_assert(is_power_2(Alignment), "Aligned allocator must be a power of 2");
public:

    typedef void value_type;
    typedef void* pointer;
    typedef const void* const_pointer;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    template <typename Other>
    struct rebind {
        typedef AlignedAllocator<Other, Alignment> other;
    };

    template <typename Other>
    constexpr AlignedAllocator& operator=(const AlignedAllocator<Other, Alignment>&) {
        return *this;
    }

};


template <size_t Alignment>
class AlignedAllocator<const void, Alignment> {
    static_assert(is_power_2(Alignment), "Aligned allocator must be a power of 2");
public:

    typedef const void value_type;
    typedef const void* pointer;
    typedef const void* const_pointer;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    template <typename Other>
    struct rebind {
        typedef AlignedAllocator<Other, Alignment> other;
    };

    template <typename Other>
    constexpr AlignedAllocator& operator=(const AlignedAllocator<Other, Alignment>&) {
        return *this;
    }

};



template <typename Sy, size_t SAlign, typename Ty, size_t TAlign>
constexpr bool operator==(const AlignedAllocator<Sy, SAlign>&, const AlignedAllocator<Ty, TAlign>&) noexcept {
    return SAlign == TAlign;
}

template <typename Sy, size_t SAlign, typename Ty, size_t TAlign>
constexpr bool operator!=(const AlignedAllocator<Sy, SAlign>&, const AlignedAllocator<Ty, TAlign>&) noexcept {
    return SAlign != TAlign;
}


}


template<typename T>
using aligned_vec = std::vector<T, dtl::AlignedAllocator<T, 16>>;

} // namespace recombine

#endif // ALIGNED_VEC_H
