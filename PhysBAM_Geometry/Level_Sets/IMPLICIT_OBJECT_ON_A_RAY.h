//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_ON_A_RAY
//##################################################################### 
#ifndef __IMPLICIT_OBJECT_ON_A_RAY__
#define __IMPLICIT_OBJECT_ON_A_RAY__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
namespace PhysBAM{

template<class TV> class RAY;

template<class T_IMPLICIT_OBJECT>
class IMPLICIT_OBJECT_ON_A_RAY:public NONLINEAR_FUNCTION<typename T_IMPLICIT_OBJECT::VECTOR_T::SCALAR(typename T_IMPLICIT_OBJECT::VECTOR_T::SCALAR)>
{
    typedef typename T_IMPLICIT_OBJECT::VECTOR_T TV;
    typedef typename TV::SCALAR T;
public:
    const T_IMPLICIT_OBJECT& implicit_object;
    RAY<TV>& ray;
        
    IMPLICIT_OBJECT_ON_A_RAY(const T_IMPLICIT_OBJECT& implicit_object_input,RAY<TV>& ray_input)
        :implicit_object(implicit_object_input),ray(ray_input)
    {}

    T operator()(const T x) const PHYSBAM_OVERRIDE
    {return implicit_object(ray.Point(x));}

//#####################################################################
};   
}
#endif


