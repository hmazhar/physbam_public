//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef _COMPONENTWISE_MIN_MAX_
#define _COMPONENTWISE_MIN_MAX_
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
namespace PhysBAM{
//#####################################################################
// Function Componentwise_Min
//#####################################################################
template<class T>
inline T Componentwise_Min(const T& v1,const T& v2)
{return T::Componentwise_Min(v1,v2);}

inline float Componentwise_Min(const float& v1,const float& v2)
{return min(v1,v2);}

inline double Componentwise_Min(const double& v1,const double& v2)
{return min(v1,v2);}

template<class T>
inline T Componentwise_Min(const T& v1,const T& v2,const T& v3)
{return Componentwise_Min(v1,Componentwise_Min(v2,v3));}

template<class T>
inline T Componentwise_Min(const T& v1,const T& v2,const T& v3,const T& v4)
{return Componentwise_Min(v1,Componentwise_Min(v2,v3,v4));}

template<class T>
inline T Componentwise_Min(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5)
{return Componentwise_Min(v1,Componentwise_Min(v2,v3,v4,v5));}

template<class T>
inline T Componentwise_Min(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6)
{return Componentwise_Min(v1,Componentwise_Min(v2,v3,v4,v5,v6));}

template<class T>
inline T Componentwise_Min(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7)
{return Componentwise_Min(v1,Componentwise_Min(v2,v3,v4,v5,v6,v7));}

template<class T>
inline T Componentwise_Min(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7,const T& v8)
{return Componentwise_Min(v1,Componentwise_Min(v2,v3,v4,v5,v6,v7,v8));}

//#####################################################################
// Function Componentwise_Max
//#####################################################################
template<class T>
inline T Componentwise_Max(const T& v1,const T& v2)
{return T::Componentwise_Max(v1,v2);}

inline float Componentwise_Max(const float& v1,const float& v2)
{return max(v1,v2);}

inline double Componentwise_Max(const double& v1,const double& v2)
{return max(v1,v2);}

template<class T>
inline T Componentwise_Max(const T& v1,const T& v2,const T& v3)
{return Componentwise_Max(v1,Componentwise_Max(v2,v3));}

template<class T>
inline T Componentwise_Max(const T& v1,const T& v2,const T& v3,const T& v4)
{return Componentwise_Max(v1,Componentwise_Max(v2,v3,v4));}

template<class T>
inline T Componentwise_Max(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5)
{return Componentwise_Max(v1,Componentwise_Max(v2,v3,v4,v5));}

template<class T>
inline T Componentwise_Max(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6)
{return Componentwise_Max(v1,Componentwise_Max(v2,v3,v4,v5,v6));}

template<class T>
inline T Componentwise_Max(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7)
{return Componentwise_Max(v1,Componentwise_Max(v2,v3,v4,v5,v6,v7));}

template<class T>
inline T Componentwise_Max(const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7,const T& v8)
{return Componentwise_Max(v1,Componentwise_Max(v2,v3,v4,v5,v6,v7,v8));}

}
#endif
