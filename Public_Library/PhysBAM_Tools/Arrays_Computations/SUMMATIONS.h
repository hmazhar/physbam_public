//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SUMMATIONS
//#####################################################################
#ifndef __SUMMATIONS__
#define __SUMMATIONS__

#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <cassert>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID>  class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class ID>
    T Sum(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result+=self(i);return result;}

    template<class T,class T_ARRAY,class ID>
    double Sum_Double_Precision(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();double result=0;ID m=self.Size();for(ID i(1);i<=m;i++) result+=self(i);return result;}

    template<class T,class T_ARRAY,class ID>
    T Sumabs(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result+=abs(self(i));return result;}
    
    template<class T,class T_ARRAY,class ID>
    T Average(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();return self.Size()?Sum(a)/typename ARRAY_BASE<T,T_ARRAY,ID>::SCALAR(self.Size()):T();}
    
    template<class T_ARRAY1,class T_ARRAY2>
    typename T_ARRAY2::ELEMENT Weighted_Sum(const T_ARRAY1& weights,const T_ARRAY2& array)
    {STATIC_ASSERT_SAME(typename T_ARRAY1::ELEMENT,typename T_ARRAY2::SCALAR);assert(weights.Size()==array.Size());
    typename T_ARRAY2::ELEMENT result((typename T_ARRAY2::ELEMENT()));typename T_ARRAY2::INDEX m=array.Size();for(typename T_ARRAY2::INDEX i(1);i<=m;i++) result+=weights(i)*array(i);return result;}
}
}
#endif
