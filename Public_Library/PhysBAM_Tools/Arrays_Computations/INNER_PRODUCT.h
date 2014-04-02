//#####################################################################
// Copyright 2009, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INNER_PRODUCT
//#####################################################################
#ifndef __INNER_PRODUCT__
#define __INNER_PRODUCT__

#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class T_ARRAY2,class ID> typename T_ARRAY2::SCALAR
    Inner_Product(const ARRAY_BASE<typename T_ARRAY2::SCALAR,T_ARRAY,ID>& m,const ARRAY_BASE<T,T_ARRAY2,ID>& a1,const ARRAY_BASE<T,T_ARRAY2,ID>& a2)
    {assert(a1.Size()==a2.Size());return ARRAYS_COMPUTATIONS::Sum((m*(a1*a2)));}

    template<class T,class T2,class T_ARRAY,class T_ARRAY2,class ID> typename T_ARRAY2::ELEMENT::SCALAR
    Inner_Product(const ARRAY_BASE<T,T_ARRAY,ID>& m,const ARRAY_BASE<T2,T_ARRAY2,ID>& a1,const ARRAY_BASE<T2,T_ARRAY2,ID>& a2)
    {assert(a1.Size()==a2.Size());
    typename T_ARRAY2::SCALAR result(0);ID size=a1.Size();for(ID i(1);i<=size;i++) result+=m(i).Inner_Product(a1(i),a2(i));return result;}

    template<class T,class T_ARRAY,class T_ARRAY2,class ID> double
    Inner_Product_Double_Precision(const ARRAY_BASE<typename T_ARRAY2::SCALAR,T_ARRAY,ID>& m,const ARRAY_BASE<T,T_ARRAY2,ID>& a1,const ARRAY_BASE<T,T_ARRAY2,ID>& a2)
    {assert(a1.Size()==a2.Size());return ARRAYS_COMPUTATIONS::Sum((m*(a1*a2))).Sum();}

template<class T,class T2,class T_ARRAY,class T_ARRAY2,class ID> typename DISABLE_IF<IS_SCALAR<T>::value,double>::TYPE
    Inner_Product_Double_Precision(const ARRAY_BASE<T,T_ARRAY,ID>& m,const ARRAY_BASE<T2,T_ARRAY2,ID>& a1,const ARRAY_BASE<T2,T_ARRAY2,ID>& a2)
    {assert(a1.Size()==a2.Size());
    double result(0);ID size=a1.Size();for(ID i(1);i<=size;i++) result+=m(i).Inner_Product(a1(i),a2(i));return result;}
}
}
#endif
