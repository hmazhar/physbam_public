//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE__
#define __LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE__

#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_RLE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_RLE<T_GRID>
class LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE:public INTERPOLATION_RLE<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:
    LINEAR_INTERPOLATION_RLE<T_GRID,T2,T_FACE_LOOKUP> interpolation;
    
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE()
    {}

    TV From_Block_Face(const T_BLOCK& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(block.Inside(X));
    u.Set_Reference_Point(block,X);TV result=T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,u,X);u.Clear_Reference_Point();
    return result;}

    T From_Block_Face_Component(const int axis,const T_BLOCK& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE
    {assert(block.Inside(X));
    u.Set_Reference_Point(block,X);T result=T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Component(axis,block,u,X);u.Clear_Reference_Point();
    return result;}

protected:
    virtual T2 Invalid_Value_Replacement(const ARRAY<T2>& u,const TV& X,const T_BLOCK& block,const TV& intersection_point,
        const int body_id,const int aggregate_id,bool& valid,const T ray_t_max=0) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

//#####################################################################
};
}
#endif
#endif
