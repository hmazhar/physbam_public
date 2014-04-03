//#####################################################################
// Copyright 2002-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_LEVELSET_ON_A_RAY
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_LEVELSET_ON_A_RAY__
#define __RLE_LEVELSET_ON_A_RAY__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class RLE_LEVELSET_ON_A_RAY:public NONLINEAR_FUNCTION<typename TV::SCALAR(typename TV::SCALAR)>
{
    typedef typename TV::SCALAR T;typedef typename RLE_GRID_POLICY<TV>::RLE_GRID T_GRID;typedef typename T_GRID::BLOCK T_BLOCK;
public:
    const LEVELSET_RLE<T_GRID>& levelset;
    RAY<TV>& ray;
    mutable T_BLOCK* last_block;

    RLE_LEVELSET_ON_A_RAY(const LEVELSET_RLE<T_GRID>& levelset_input,RAY<TV>& ray_input)
        :levelset(levelset_input),ray(ray_input),last_block(0)
    {}

    RLE_LEVELSET_ON_A_RAY(const RLE_IMPLICIT_OBJECT<TV>& levelset_input,RAY<TV>& ray_input)
        :levelset(levelset_input.levelset),ray(ray_input),last_block(0)
    {}

    ~RLE_LEVELSET_ON_A_RAY()
    {delete last_block;}

    T operator()(const T x) const PHYSBAM_OVERRIDE
    {TV X=ray.Point(x);
    if(!last_block) last_block=new T_BLOCK(levelset.grid,X);
    else if(!*last_block || !last_block->Lazy_Inside(X)) last_block->Initialize(X);
    return levelset.Phi_Far(*last_block,X);}

//#####################################################################
};
}
#endif
#endif
