//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CENTER
//#####################################################################
#ifndef __CENTER__
#define __CENTER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
namespace PhysBAM{

namespace POINT_CLOUDS_COMPUTATIONS
{
    template<class TV,class T_ARRAY>
    static TV Center(const ARRAY_BASE<TV,T_ARRAY>& X)
    {typename TV::SCALAR total=X.Size();return total?ARRAYS_COMPUTATIONS::Sum(X)/total:TV();}

    template<class TV,class T_ARRAY,class T_ARRAY2>
    static TV Weighted_Center(const ARRAY_BASE<TV,T_ARRAY>& X,const ARRAY_BASE<typename TV::SCALAR,T_ARRAY2>& weights)
    {typename TV::SCALAR total=ARRAYS_COMPUTATIONS::Sum(weights);return total?ARRAYS_COMPUTATIONS::Weighted_Sum(weights,X)/total:TV();}
}
}
#endif
