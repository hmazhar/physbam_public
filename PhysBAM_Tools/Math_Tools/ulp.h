//#####################################################################
// Copyright 2008, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function ulp
//#####################################################################
#ifndef __ulp__
#define __ulp__
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <cassert>
#include <limits>

namespace PhysBAM{

inline float ulp(float x)
{STATIC_ASSERT((sizeof(float)==sizeof(unsigned int)&&std::numeric_limits<float>::is_iec559));
assert(x!=0);
union {unsigned int i;float f;} u;
u.f=x;
return ((u.i^=1)&1?u.f-x:x-u.f);}

inline double ulp(double x)
{STATIC_ASSERT((sizeof(double)==sizeof(unsigned long long)&&std::numeric_limits<double>::is_iec559));
assert(x!=0);union {unsigned long long i;double f;} u;u.f=x;
return ((u.i^=1)&1?u.f-x:x-u.f);}

}
#endif
