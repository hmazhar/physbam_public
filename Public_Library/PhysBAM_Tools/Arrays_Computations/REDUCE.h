//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REDUCE
//#####################################################################
#ifndef __REDUCE__
#define __REDUCE__

namespace PhysBAM{

namespace ARRAYS_COMPUTATIONS
{
    template<class T_ARRAY> static typename T_ARRAY::ELEMENT Reduce(typename T_ARRAY::ELEMENT (*reduce_func)(typename T_ARRAY::ELEMENT,typename T_ARRAY::ELEMENT),const T_ARRAY& array)
    {typedef typename T_ARRAY::ELEMENT T;
    T sum=T();for(int i=1;i<=array.Size();i++) sum=reduce_func(sum,array(i));
    return sum;}
}
}
#endif
