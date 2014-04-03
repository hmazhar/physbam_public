//#####################################################################
// Copyright 2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function rint
//#####################################################################
#ifndef __rint__
#define __rint__

#include <cmath>
namespace PhysBAM{

#ifdef WIN32
inline float rint(const float x){return floorf(x+(x>0?.5f:-.5f));}
inline double rint(const double x){return floor(x+(x>0?.5:-.5));}
#else
inline float rint(const float x){return ::rintf(x);}
inline double rint(const double x){return ::rint(x);}
#endif

}
#endif
