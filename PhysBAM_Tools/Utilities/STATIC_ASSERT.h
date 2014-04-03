//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header STATIC_ASSERT
//#####################################################################
#ifndef __STATIC_ASSERT__
#define __STATIC_ASSERT__

namespace PhysBAM{
#define PHYSBAM_JOIN( X, Y ) PHYSBAM_DO_JOIN( X, Y )
#define PHYSBAM_DO_JOIN( X, Y ) PHYSBAM_DO_JOIN2(X,Y)
#define PHYSBAM_DO_JOIN2( X, Y ) X##Y

template <bool x> struct STATIC_ASSERTION_FAILURE;

template <> struct STATIC_ASSERTION_FAILURE<true> { enum { value = 1 }; };

template<int x> struct static_assert_test{};

#define PHYSBAM_STATIC_ASSERT( B ) \
   typedef static_assert_test<\
      sizeof(STATIC_ASSERTION_FAILURE< (bool)( B ) >)>\
         PHYSBAM_JOIN(boost_static_assert_typedef_, __LINE__)

#define STATIC_ASSERT(...) PHYSBAM_STATIC_ASSERT((__VA_ARGS__))
}

#endif
