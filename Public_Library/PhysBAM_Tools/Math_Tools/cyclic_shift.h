//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function cyclic shifts
//#####################################################################
//
// Cyclic Shifts 3 or 4 values
//               
//#####################################################################
#ifndef __cyclic_shift__
#define __cyclic_shift__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> inline void
cyclic_shift(T& i,T& j,T& k) 
{T temp=k;k=j;j=i;i=temp;}

template<class T> inline void
cyclic_shift(T& i,T& j,T& k,T& l) 
{T temp=l;l=k;k=j;j=i;i=temp;}

template<class T> inline void
cyclic_shift(VECTOR<T,1>& v)
{}

template<class T> inline void
cyclic_shift(VECTOR<T,2>& v)
{exchange(v.x,v.y);}

template<class T> inline void
cyclic_shift(VECTOR<T,3>& v)
{cyclic_shift(v.x,v.y,v.z);}

template<class T> inline void
cyclic_shift(VECTOR<T,4>& v)
{cyclic_shift(v[1],v[2],v[3],v[4]);}

}
#endif
