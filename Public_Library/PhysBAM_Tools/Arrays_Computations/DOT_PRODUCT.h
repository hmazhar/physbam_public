//#####################################################################
// Copyright 2009, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DOT_PRODUCT
//#####################################################################
#ifndef __DOT_PRODUCT__
#define __DOT_PRODUCT__

#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <cassert>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class T_ARRAY2,class ID> typename SCALAR_POLICY<T>::TYPE
    Dot_Product(const ARRAY_BASE<T,T_ARRAY,ID>& a1,const ARRAY_BASE<T,T_ARRAY2,ID>& a2)
    {assert(a1.Size()==a2.Size());
    typename SCALAR_POLICY<T>::TYPE result(0);ID m=a1.Size();for(ID i(1);i<=m;i++) result+=PhysBAM::Dot_Product(a1(i),a2(i));return result;}

    template<class T,class T_ARRAY,class T_ARRAY2,class ID> double
    Dot_Product_Double_Precision(const ARRAY_BASE<T,T_ARRAY,ID>& a1,const ARRAY_BASE<T,T_ARRAY2,ID>& a2)
    {assert(a1.Size()==a2.Size());
    double result(0);ID m=a1.Size();for(ID i(1);i<=m;i++) result+=PhysBAM::Dot_Product(a1(i),a2(i));return result;}
}
}
#endif
