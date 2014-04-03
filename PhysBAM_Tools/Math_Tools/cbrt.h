//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function cbrt
//#####################################################################
//
// Apparently cbrt is standard C but not standard C++.
//               
//#####################################################################
#ifndef __cbrt__
#define __cbrt__

#include <cmath>
namespace PhysBAM{

#ifdef WIN32
inline float cbrt(const float a){return std::pow(a,(float)(1/3.));}
inline double cbrt(const double a){return std::pow(a,(double)(1/3.));}
#else
inline float cbrt(const float a){return ::cbrtf(a);}
inline double cbrt(const double a){return ::cbrt(a);}
#endif

}
#endif
