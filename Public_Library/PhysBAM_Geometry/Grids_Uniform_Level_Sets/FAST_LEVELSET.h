//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_LEVELSET
//#####################################################################
#ifndef __FAST_LEVELSET__
#define __FAST_LEVELSET__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID_input>
class FAST_LEVELSET:public LEVELSET_POLICY<T_GRID_input>::LEVELSET
{
    typedef T_GRID_input T_GRID;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::AVERAGING T_AVERAGING;
public:
    typedef T_LEVELSET BASE;
    using BASE::grid;using BASE::boundary;using BASE::phi;using BASE::curvature_motion;using BASE::sigma;using BASE::max_time_step;
    using BASE::small_number;using BASE::Set_Face_Velocities_Valid_Mask;

    T half_band_width;

    FAST_LEVELSET(T_GRID& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells=3);
    ~FAST_LEVELSET();

    void Set_Band_Width(const T number_of_cells=6)
    {half_band_width=number_of_cells*grid.dX.Max()/2;}

//#####################################################################
    T CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const;
    T CFL(const T_ARRAYS_VECTOR& velocity) const;
public:
    void Fast_Marching_Method(const int local_advection_spatial_order, const T time=0);
    void Fast_Marching_Method_Outside_Band(const int local_advection_spatial_order, const T time=0);
//#####################################################################
};
}
#endif
