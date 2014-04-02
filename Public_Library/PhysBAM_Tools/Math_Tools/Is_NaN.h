//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class Is_NaN
//#####################################################################
#ifndef __Is_NaN__
#define __Is_NaN__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T>
inline bool Is_NaN(const T& scalar)
{
#ifdef WIN32
    return _isnan(scalar);
#elif defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)
    return isnan(scalar);
#else
    return false;
#endif
}

template<class T>
inline bool Is_NaN(const VECTOR<T,1>& v)
{return Is_NaN(v.x);}

template<class T>
inline bool Is_NaN(const VECTOR<T,2>& v)
{return Is_NaN(v.x)||Is_NaN(v.y);}

template<class T>
inline bool Is_NaN(const VECTOR<T,3>& v)
{return Is_NaN(v.x)||Is_NaN(v.y)||Is_NaN(v.z);}

}
#endif
