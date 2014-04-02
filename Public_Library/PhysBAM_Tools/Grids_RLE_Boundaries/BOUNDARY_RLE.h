//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY
//#####################################################################
//
// Inherited by all BOUNDARY classes.
// The ghost cell routines copy u into u_ghost and add a layer of three ghost cells around the boundary of u.
// The boundary condition routines modify u directly.
//
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __BOUNDARY_RLE__
#define __BOUNDARY_RLE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_POLICY_RLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T_GRID,class T_ITERATOR> class RLE_GRID_TRANSFER_ITERATOR;

template<class T_GRID,class T2>
class BOUNDARY_RLE:public BOUNDARY<typename T_GRID::VECTOR_T,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2_BASE;
    typedef typename T_FACE_ARRAYS::template REBIND<T2>::TYPE T_FACE_ARRAYS_T2_BASE;
    typedef typename REBIND<typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR,T2>::TYPE T_BOUNDARY_T2;
public:
    typedef BOUNDARY<typename T_GRID::VECTOR_T,T2> BASE;
    using BASE::use_fixed_boundary;using BASE::clamp_below;using BASE::clamp_above;using BASE::fixed_boundary_value;using BASE::lower_threshold;using BASE::upper_threshold;
    using BASE::Set_Constant_Extrapolation;

    BOUNDARY_RLE();
    BOUNDARY_RLE(const TV_SIDES& constant_extrapolation);
    virtual ~BOUNDARY_RLE();

    void Fill_Ghost_Cells_Cell(const T_GRID& grid,const T_ARRAYS_T2_BASE& u,T_ARRAYS_T2_BASE& u_ghost,const T time)
    {Fill_Ghost_Cells(grid,u,u_ghost,0,time);}

//#####################################################################
    virtual void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2_BASE& u,const T time){}
    virtual void Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2_BASE& u,const T time){}
    virtual void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2_BASE& u,T_ARRAYS_T2_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3);
    virtual void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2_BASE& u,T_FACE_ARRAYS_T2_BASE& u_ghost,const T time,const int number_of_ghost_cells=3);
    template<class T3,class T4,class T_ITERATOR> void
    Fill_Ghost_Cells_Helper(const T_GRID& grid,ARRAY_BASE<T4,ARRAY<T4> >& u_ghost,const T3& fixed_boundary_value,const T3& lower_threshold,
        const T3& upper_threshold,const int number_of_ghost_cells) const;
    template<class T3,class T4,class T_ITERATOR> void Fill_Single_Ghost_Region(const T_GRID& grid,ARRAY_BASE<T4,ARRAY<T4> >& u_ghost,const int horizontal_side,
        const typename T_GRID::BOX_HORIZONTAL_INT& region,const T3& fixed_boundary_value,const T3& lower_threshold,
        const typename ENABLE_IF<IS_SAME<T3,T4>::value,T3>::TYPE& upper_threshold) const;
    template<class T3,class T4,class T_ITERATOR> void Fill_Single_Ghost_Region(const T_GRID& grid,ARRAY_BASE<T4,ARRAY<T4> >& u_ghost,const int horizontal_side,
        const typename T_GRID::BOX_HORIZONTAL_INT& region,const T3& fixed_boundary_value,const T3& lower_threshold,
        const typename DISABLE_IF<IS_SAME<T3,T4>::value,T3>::TYPE& upper_threshold) const;
//#####################################################################
};
}
#endif
#endif
