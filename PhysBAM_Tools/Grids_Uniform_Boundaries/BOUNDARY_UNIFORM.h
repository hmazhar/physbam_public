//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_UNIFORM
//#####################################################################
//
// The ghost cell routines copy u into u_ghost and add a layer of three ghost cells around the boundary of u.
// The boundary condition routines modify u directly.
//
//#####################################################################
#ifndef __BOUNDARY_UNIFORM__
#define __BOUNDARY_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_UNIFORM:public BOUNDARY<typename T_GRID::VECTOR_T,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,T2>::TYPE T_FACE_ARRAYS_T2;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    typedef BOUNDARY<typename T_GRID::VECTOR_T,T2> BASE;
    using BASE::use_fixed_boundary;using BASE::clamp_below;using BASE::clamp_above;using BASE::fixed_boundary_value;using BASE::lower_threshold;using BASE::upper_threshold;

    BOUNDARY_UNIFORM();
    virtual ~BOUNDARY_UNIFORM();

protected:
    int Boundary(const int side,const RANGE<TV_INT>& region) const
    {int axis=(side+1)/2;return side&1?region.Maximum_Corner()[axis]+1:region.Minimum_Corner()[axis]-1;}
public:

    void Fill_Ghost_Cells_Cell(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells)
    {Fill_Ghost_Cells(grid,u,u_ghost,0,time,number_of_ghost_cells);}

//#####################################################################
    virtual void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time);
    virtual void Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2& u,const T time);
    virtual void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells);
    virtual void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells);
    virtual void Apply_Boundary_Condition_Single_Side(const T_GRID& grid,T_ARRAYS_T2& u,const int side,const T time) const;
    virtual void Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const;
    void Find_Ghost_Regions(const T_GRID& grid,ARRAY<RANGE<TV_INT> >& regions,const int number_of_ghost_cells) const;
    void Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const int side,const RANGE<TV_INT>& region) const;
    void Fill_Single_Ghost_Region_Threaded(RANGE<TV_INT>& region,const T_GRID& grid,T_ARRAYS_T2& u_ghost,const int side);
//#####################################################################
};
}
#endif
