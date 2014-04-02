//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LINEAR_INTERPOLATION_RLE__
#define __LINEAR_INTERPOLATION_RLE__

#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_RLE<T_GRID>
class LINEAR_INTERPOLATION_RLE:public INTERPOLATION_RLE<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename T_GRID::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:
    template<class T3> struct REBIND{typedef LINEAR_INTERPOLATION_RLE<T_GRID,T3,T_FACE_LOOKUP> TYPE;};

    T2 From_Block_Cell(const T_BLOCK& block,const ARRAY<T2>& u,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Cell(block,u,X);}

    TV From_Block_Face(const T_BLOCK& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,u,X);}

    T From_Block_Face_Component(const int axis,const T_BLOCK& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Component(axis,block,u,X);}

//#####################################################################
};
}
#endif
#endif
