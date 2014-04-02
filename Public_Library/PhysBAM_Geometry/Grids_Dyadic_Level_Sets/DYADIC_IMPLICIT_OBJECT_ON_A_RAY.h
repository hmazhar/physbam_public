//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_IMPLICIT_SURFACE_ON_A_RAY
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __DYADIC_IMPLICIT_OBJECT_ON_A_RAY__
#define __DYADIC_IMPLICIT_OBJECT_ON_A_RAY__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class T_GRID>
class DYADIC_IMPLICIT_OBJECT_ON_A_RAY:public NONLINEAR_FUNCTION<typename T_GRID::SCALAR(typename T_GRID::SCALAR)>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::CELL T_CELL;
public:
    const DYADIC_IMPLICIT_OBJECT<TV>& implicit_object;
    RAY<TV>& ray;
    mutable BLOCK_DYADIC<T_GRID>* last_block;
    mutable const ARRAY<T>* phi_nodes;

    DYADIC_IMPLICIT_OBJECT_ON_A_RAY(const DYADIC_IMPLICIT_OBJECT<TV>& implicit_object_input,RAY<TV>& ray_input)
        :implicit_object(implicit_object_input),ray(ray_input),last_block(0),phi_nodes(0)
    {}

    ~DYADIC_IMPLICIT_OBJECT_ON_A_RAY()
    {delete last_block;}

    T operator()(const T x) const PHYSBAM_OVERRIDE
    {TV X=ray.Point(x);
    const T_CELL* base_cell=0;
    if(last_block)base_cell=implicit_object.levelset.grid.Base_Cell_By_Neighbor_Path(last_block->Base_Cell(),X);
    else base_cell=implicit_object.levelset.grid.Base_Cell(X);
    if(base_cell){
        if(!last_block) last_block=new BLOCK_DYADIC<T_GRID>(implicit_object.levelset.grid,base_cell);
        else if(base_cell!=last_block->Base_Cell()) last_block->Initialize(base_cell);
        return implicit_object.levelset.Phi(*last_block,X);}
    else{assert(phi_nodes);return implicit_object.levelset.Phi(X,phi_nodes);}} // if it's called outside the band, we may have to fall back on nodal interpolation

//#####################################################################
};
}
#endif
#endif
