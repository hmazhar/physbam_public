//#####################################################################
// Copyright 2005, Jiayi Chong, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MULTIPLE_ON_A_RAY
//#####################################################################
#ifndef __LEVELSET_MULTIPLE_ON_A_RAY__
#define __LEVELSET_MULTIPLE_ON_A_RAY__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE.h>
namespace PhysBAM{

template<class T_LEVELSET_MULTIPLE>
class LEVELSET_MULTIPLE_ON_A_RAY:public NONLINEAR_FUNCTION<typename T_LEVELSET_MULTIPLE::VECTOR_T::SCALAR(typename T_LEVELSET_MULTIPLE::VECTOR_T::SCALAR)>
{
    typedef typename T_LEVELSET_MULTIPLE::VECTOR_T TV;
    typedef typename TV::SCALAR T;
public:
    const T_LEVELSET_MULTIPLE& levelset_multiple;
    RAY<VECTOR<T,3> >& ray;
    int region;

    LEVELSET_MULTIPLE_ON_A_RAY(const T_LEVELSET_MULTIPLE& levelset_multiple_input,RAY<VECTOR<T,3> >& ray_input,const int region_input)
        :levelset_multiple(levelset_multiple_input),ray(ray_input),region(region_input)
    {}

    T operator()(const T x) const PHYSBAM_OVERRIDE
    {VECTOR<T,3> X=ray.Point(x);int minimum_region,second_minimum_region;T minimum_phi,second_minimum_phi;
    levelset_multiple.Two_Minimum_Regions(X,minimum_region,second_minimum_region,minimum_phi,second_minimum_phi);
    return levelset_multiple.Phi(region,X)-(T).5*(minimum_phi+second_minimum_phi);}

//#####################################################################
};
}
#endif
