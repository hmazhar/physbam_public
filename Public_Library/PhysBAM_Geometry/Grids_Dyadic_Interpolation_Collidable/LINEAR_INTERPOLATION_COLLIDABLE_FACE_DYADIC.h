//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC__
#define __LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC__

#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID>
class LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC:public INTERPOLATION_DYADIC<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL T_CELL;typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:
    LINEAR_INTERPOLATION_DYADIC<T_GRID,T2,T_FACE_LOOKUP> interpolation;
    
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC()
    {}

    TV From_Block_Face(const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(block.Inside(X));
    u.Set_Reference_Point(X);TV result=T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,u,X);u.Clear_Reference_Point();
    return result;}

    T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE
    {assert(block.Inside(X));
    u.Set_Reference_Point(X);T result=T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Component(axis+1,block,u,X);u.Clear_Reference_Point();
    return result;}

protected:
    virtual T2 Invalid_Value_Replacement(const T_GRID& grid,const ARRAY<T2>& u,const TV& X,const BLOCK_DYADIC<T_GRID>& block,const TV& intersection_point,
        const int body_id,const int aggregate_id,bool& valid,const T ray_t_max=0) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

//#####################################################################
};
}
#endif
#endif
