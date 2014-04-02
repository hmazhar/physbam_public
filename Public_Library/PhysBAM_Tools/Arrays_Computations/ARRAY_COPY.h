//#####################################################################
// Copyright 2009, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COPY
//#####################################################################
#ifndef __ARRAY_COPY__
#define __ARRAY_COPY__

#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <cassert>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID>  class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class ID>
    void Fill(ARRAY_BASE<T,T_ARRAY,ID>& a,const T& constant)
    {T_ARRAY& self=a.Derived();ID m=self.Size();for(ID i(1);i<=m;i++) self(i)=constant;}

    template<class T,class T2,class T_ARRAY1,class T_ARRAY2,class ID>
    void Copy_With_Offset(const ARRAY_BASE<T,T_ARRAY1,ID>& old_copy,ARRAY_BASE<T2,T_ARRAY2,ID>& new_copy,const ID offset)
    {STATIC_ASSERT(CAN_ASSIGN<T,T2>::value);
    ID m=old_copy.Size();assert(m+offset<=new_copy.Size());
    for(ID i(1);i<=m;i++) new_copy(i+offset)=old_copy(i);}
}
}
#endif
