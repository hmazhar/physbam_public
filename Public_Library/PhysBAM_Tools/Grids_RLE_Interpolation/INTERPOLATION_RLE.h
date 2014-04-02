//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERPOLATION_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __INTERPOLATION_RLE__
#define __INTERPOLATION_RLE__

#include <PhysBAM_Tools/Grids_RLE_Interpolation/FACE_LOOKUP_RLE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_RLE<T_GRID>
class INTERPOLATION_RLE:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::BLOCK T_BLOCK;
public:
    template<class T3> struct REBIND{typedef INTERPOLATION_RLE<T_GRID,T3,T_FACE_LOOKUP> TYPE;};

    virtual ~INTERPOLATION_RLE()
    {}

//#####################################################################
    virtual T2 From_Block_Cell(const T_BLOCK& block,const ARRAY<T2>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV From_Block_Face(const T_BLOCK& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T From_Block_Face_Component(const int axis,const T_BLOCK& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
#endif
