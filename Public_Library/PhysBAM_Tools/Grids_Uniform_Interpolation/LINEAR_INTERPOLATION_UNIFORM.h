//#####################################################################
// Copyright 2005-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_UNIFORM__
#define __LINEAR_INTERPOLATION_UNIFORM__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_1D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/Componentwise_Min_Max.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID>
class LINEAR_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> > T_ARRAYS_T2;
    //typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
public:
    template<class T3> struct REBIND{typedef LINEAR_INTERPOLATION_UNIFORM<T_GRID,T3,T_FACE_LOOKUP> TYPE;};
    
    T2 Clamped_To_Array_No_Extrema(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
    // requires function to be monotonic increasing; clamps values out of range
    static void Compute_Inverse_Map(const GRID<TV>& domain_grid,const ARRAY<T2,VECTOR<int,1> >& function,const GRID<TV>& range_grid,ARRAY<T2,VECTOR<int,1> >& inverse_function);
    VECTOR<T2,2> Extrema_Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const  PHYSBAM_OVERRIDE;
    VECTOR<T2,2> Extrema_From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;

    TV From_Block_Face(const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,u,X);}

    T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > From_Block_Face_Component_Weights(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE;

    VECTOR<TV,2> Extrema_From_Block_Face(const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Extrema_Face(block,u_min,u_max,X);}

    VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Extrema_Face_Component(axis,block,u_min,u_max,X);}

//#####################################################################
};
}
#endif
