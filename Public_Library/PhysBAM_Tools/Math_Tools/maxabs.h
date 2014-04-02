//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function maxabs
//#####################################################################
//
// finds the maximum absolute value
//
//#####################################################################
#ifndef __maxabs__
#define __maxabs__

#include <PhysBAM_Tools/Math_Tools/max.h>
#include <cmath>
namespace PhysBAM{

using std::abs;

template<class T>
inline T maxabs(const T a,const T b)
{return max(abs(a),abs(b));}

template<class T>
inline T maxabs(const T a,const T b,const T c)
{return max(maxabs(a,b),abs(c));}

template<class T>
inline T maxabs(const T a,const T b,const T c,const T d)
{return max(maxabs(a,b,c),abs(d));}

template<class T>
inline T maxabs(const T a,const T b,const T c,const T d,const T e)
{return max(maxabs(a,b,c,d),abs(e));}

template<class T>
inline T maxabs(const T a,const T b,const T c,const T d,const T e,const T f)
{return max(maxabs(a,b,c,d,e),abs(f));}

template<class T>
inline T maxabs(const T a,const T b,const T c,const T d,const T e,const T f,const T g)
{return max(maxabs(a,b,c,d,e,f),abs(g));}

template<class T>
inline T maxabs(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h)
{return max(maxabs(a,b,c,d,e,f,g),abs(h));}

template<class T>
inline T maxabs(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h,const T i)
{return max(maxabs(a,b,c,d,e,f,g,h),abs(i));}

}
#endif
