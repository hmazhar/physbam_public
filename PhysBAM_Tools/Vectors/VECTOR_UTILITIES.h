//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VECTOR_UTILITIES__
#define __VECTOR_UTILITIES__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

namespace VECTOR_UTILITIES
{
//#####################################################################
// Function Complement
//#####################################################################
inline bool Complement(const bool b)
{return !b;}

template<class T,int d> VECTOR<T,d> Complement(const VECTOR<T,d>& v)
{VECTOR<T,d> v_complement;for(int i=1;i<=d;i++) v_complement(i)=Complement(v(i));return v_complement;}

template<class T> VECTOR<T,1> Complement(const VECTOR<T,1>& v)
{return VECTOR<T,1>(Complement(v.x));}

template<class T> VECTOR<T,2> Complement(const VECTOR<T,2>& v)
{return VECTOR<T,2>(Complement(v.x),Complement(v.y));}

template<class T> VECTOR<T,3> Complement(const VECTOR<T,3>& v)
{return VECTOR<T,3>(Complement(v.x),Complement(v.y),Complement(v.z));}

}
}
#endif
