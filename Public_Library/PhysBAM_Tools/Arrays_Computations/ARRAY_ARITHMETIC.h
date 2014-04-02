//#####################################################################
// Copyright 2009, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_ARITHMETIC
//#####################################################################
#ifndef __ARRAY_ARITHMETIC__
#define __ARRAY_ARITHMETIC__

#include <PhysBAM_Tools/Math_Tools/Inverse.h>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID>  class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class ID,class T_ARRAY1> T_ARRAY&
    Plus_Equals(ARRAY_BASE<T,T_ARRAY,ID>& a,const ARRAY_BASE<T,T_ARRAY1,ID>& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();const T_ARRAY1& b_=b.Derived();assert(m==b_.Size());
    for(ID i(1);i<=m;i++) self(i)+=b_(i);return self;}

    template<class T,class T_ARRAY,class ID> T_ARRAY&
    Plus_Equals(ARRAY_BASE<T,T_ARRAY,ID>& a,const T& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)+=b;return self;}

    template<class T,class T_ARRAY,class ID,class T_ARRAY1> T_ARRAY&
    Minus_Equals(ARRAY_BASE<T,T_ARRAY,ID>&a,const ARRAY_BASE<T,T_ARRAY1,ID>& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();const T_ARRAY1& b_=b.Derived();assert(m==b_.Size());
    for(ID i(1);i<=m;i++) self(i)-=b_(i);return self;}

    template<class T,class T_ARRAY,class ID> T_ARRAY&
    Minus_Equals(ARRAY_BASE<T,T_ARRAY,ID>&a,const T& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)-=b;return self;}

    template<class T,class T_ARRAY,class ID,class T2,class T_ARRAY_T2> T_ARRAY&
    Times_Equals(ARRAY_BASE<T,T_ARRAY,ID>& a,const ARRAY_BASE<T2,T_ARRAY_T2,ID>& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();const T_ARRAY_T2& b_=b.Derived();assert(m==b_.Size());
    for(ID i(1);i<=m;i++) self(i)*=b_(i);return self;}

    template<class T,class T_ARRAY,class ID> T_ARRAY&
    Times_Equals(ARRAY_BASE<T,T_ARRAY,ID>& a,const typename T_ARRAY::SCALAR& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)*=b;return self;}

    template<class T,class T_ARRAY,class ID,class T2,class T_ARRAY_T2> T_ARRAY&
    Divide_Equals(ARRAY_BASE<T,T_ARRAY,ID>& a,const ARRAY_BASE<T2,T_ARRAY_T2,ID>& b)
    {T_ARRAY& self=a.Derived();ID m=self.Size();const T_ARRAY_T2& b_=b.Derived();assert(m==b_.Size());
    for(ID i(1);i<=m;i++){assert(b_(i));self(i)/=b_(i);}return self;}

    template<class T,class T_ARRAY,class ID> T_ARRAY&
    Divide_Equals(ARRAY_BASE<T,T_ARRAY,ID>& a,const typename T_ARRAY::SCALAR& b)
    {return Times_Equals(a,Inverse(b));}
}
}
#endif
