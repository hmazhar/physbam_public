//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM__
#define __LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID>
class LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM:public INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
public:
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP> interpolation;

    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM()
    {}

    TV From_Block_Face(const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {u.Set_Reference_Point(X);TV result=interpolation.From_Block_Face(grid,block,u,X);u.Clear_Reference_Point();return result;}

    T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE
    {u.Set_Reference_Point(X);T result=interpolation.From_Block_Face_Component(axis,grid,block,u,X);u.Clear_Reference_Point();return result;}

    VECTOR<TV,2> Extrema_From_Block_Face(const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {u_min.Set_Reference_Point(X);u_max.Set_Reference_Point(X);VECTOR<TV,2> result=interpolation.Extrema_From_Block_Face(grid,block,u_min,u_max,X);
    u_min.Clear_Reference_Point();u_max.Clear_Reference_Point();return result;}

    VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {u_min.Set_Reference_Point(X);u_max.Set_Reference_Point(X);VECTOR<T,2> result=interpolation.Extrema_From_Block_Face_Component(axis,grid,block,u_min,u_max,X);
    u_min.Clear_Reference_Point();u_max.Clear_Reference_Point();return result;}

protected:
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT // TODO(jontg): Is this in any way actually a necessary function to have?
    virtual T2 Invalid_Value_Replacement(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const BLOCK_DYADIC<T_GRID>& block,const TV& intersection_point,
        const int body_id,const int aggregate_id,bool& valid,const T ray_t_max=0) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#endif

//#####################################################################
};
}
#endif
