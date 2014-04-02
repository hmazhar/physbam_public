//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function minabs
//#####################################################################
//
// finds the minimum absolute value
//
//#####################################################################
#ifndef __minabs__
#define __minabs__

#include <PhysBAM_Tools/Math_Tools/min.h>
#include <cmath>
namespace PhysBAM{

// a should already be nonnegative
template<class T>
inline T minabs_incremental(const T a,const T b)
{return min(a,abs(b));}

template<class T>
inline T minabs(const T a,const T b)
{return min(abs(a),abs(b));}

template<class T>
inline T minabs(const T a,const T b,const T c)
{return minabs_incremental(minabs(a,b),c);}

template<class T>
inline T minabs(const T a,const T b,const T c,const T d)
{return min(minabs(a,b),minabs(c,d));}

template<class T>
inline T minabs(const T a,const T b,const T c,const T d,const T e)
{return minabs_incremental(minabs(a,b,c,d),e);}

template<class T>
inline T minabs(const T a,const T b,const T c,const T d,const T e,const T f)
{return min(minabs(a,b,c,d),minabs(e,f));}

}
#endif
