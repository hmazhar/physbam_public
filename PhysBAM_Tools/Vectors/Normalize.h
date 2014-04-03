//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Normalize 
//#####################################################################
#ifndef __Normalize__
#define __Normalize__

namespace PhysBAM{

template<class T>
inline typename T::SCALAR Normalize(T& v)
{return v.Normalize();}

template<class T>
inline T Normalized(const T& v)
{return v.Normalized();}

inline float Normalize(float& a)
{float a_save=a;
if(a>=0){a=1;return a_save;} 
else{a=-1;return -a_save;}}

inline double Normalize(double& a)
{double a_save=a;
if(a>=0){a=1;return a_save;} 
else{a=-1;return -a_save;}}

inline float Normalized(const float a)
{return a>=0?(float)1:(float)-1;}

inline double Normalized(const double a)
{return a>=0?(double)1:(double)-1;}

}
#endif
