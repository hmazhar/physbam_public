//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_DYADIC
//#####################################################################
//
// The ghost cell routines copy u into u_ghost and add a layer of three ghost cells around the boundary of u.
// The boundary condition routines modify u directly.
//
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __BOUNDARY_DYADIC__
#define __BOUNDARY_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
//#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_DYADIC:public BOUNDARY<typename T_GRID::VECTOR_T,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename T_GRID::UNIFORM_GRID T_UNIFORM_GRID;typedef typename T_GRID::CELL T_CELL;
    // typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
public:
    typedef BOUNDARY<typename T_GRID::VECTOR_T,T2> BASE;
    using BASE::use_fixed_boundary;using BASE::clamp_below;using BASE::clamp_above;using BASE::fixed_boundary_value;using BASE::lower_threshold;using BASE::upper_threshold;
    using BASE::Set_Constant_Extrapolation;

    BOUNDARY_DYADIC();
    BOUNDARY_DYADIC(const TV_SIDES& constant_extrapolation);
    virtual ~BOUNDARY_DYADIC();

//#####################################################################
    virtual void Apply_Boundary_Condition(const T_GRID& grid,ARRAY<T2>& u,const T time){}
    virtual void Fill_Ghost_Cells_Node(const T_GRID& grid,const ARRAY<T2>& u,ARRAY<T2>& u_ghost,const T time);
    virtual void Fill_Ghost_Cells_Cell(const T_GRID& grid,const ARRAY<T2>& u,ARRAY<T2>& u_ghost,const T time);
    virtual void Fill_Ghost_Cells_Face(const T_GRID& grid,const ARRAY<T>& u,ARRAY<T>& u_ghost,const T time);
protected:
    void Fill_Individual_Side_Ghost_Cells_Node(const T_GRID& grid,const int side,ARRAY<T2>& u_ghost,const T time);
    void Fill_Individual_Side_Ghost_Cells_Cell(const T_GRID& grid,const int side,ARRAY<T2>& u_ghost,const T time);
    void Fill_Individual_Side_Ghost_Cells_Face(const T_GRID& grid,const int side,ARRAY<T>& u_ghost,const T time);
//#####################################################################
};
}
#endif
#endif
