//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once
#include <stdlib.h>
//#include <malloc.h>

#if __cplusplus >= 201103L
#define NOEXCEPT_SPEC noexcept
#else
#define NOEXCEPT_SPEC
#endif
// ALIGNMENT must be a power of 2 !
#define ALIGNMENT 32

template <typename T>
class aligned_allocator {
  public:
    typedef T*              pointer;
    typedef T const*        const_pointer;
    typedef T&              reference;
    typedef T const&        const_reference;
    typedef T               value_type;
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;

    template <typename U>
    struct rebind {
        typedef aligned_allocator<U> other;
    };

    aligned_allocator() NOEXCEPT_SPEC {
    }

    aligned_allocator(aligned_allocator const& a) NOEXCEPT_SPEC {
    }

    template <typename U>
    aligned_allocator(aligned_allocator<U> const& b) NOEXCEPT_SPEC {
    }

    pointer allocate(size_type n) {
      pointer p;
      if(posix_memalign(reinterpret_cast<void**>(&p), ALIGNMENT, n*sizeof(T) ))
          throw std::bad_alloc();
      return p;
    }

    void deallocate(pointer p, size_type n) NOEXCEPT_SPEC {
        std::free(p);
    }

    size_type max_size() const NOEXCEPT_SPEC {
        std::allocator<T> a;
        return a.max_size();
    }

#if __cplusplus >= 201103L
    template <typename C, class... Args>
    void construct(C* c, Args&&... args) {
        new ((void*)c) C(std::forward<Args>(args)...);
    }
#else
    void construct(pointer p, const_reference t) {
        new((void *)p) T(t);
    }
#endif

    template <typename C>
    void destroy(C* c) {
        c->~C();
    }

    bool operator == (aligned_allocator const& a2) const NOEXCEPT_SPEC {
        return true;
    }

    bool operator != (aligned_allocator const& a2) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U>
    bool operator == (aligned_allocator<U> const& b) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U>
    bool operator != (aligned_allocator<U> const& b) const NOEXCEPT_SPEC {
        return true;
    }
};
