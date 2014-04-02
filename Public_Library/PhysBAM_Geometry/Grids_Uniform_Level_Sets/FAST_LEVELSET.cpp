//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_LEVELSET<T_GRID>::
FAST_LEVELSET(T_GRID& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input)
    :T_LEVELSET(grid_input,phi_input,number_of_ghost_cells_input)
{
    Set_Band_Width();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FAST_LEVELSET<T_GRID>::
~FAST_LEVELSET()
{}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR FAST_LEVELSET<T_GRID>::
CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    T dt_convection=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi(cell))<=half_band_width){
        T local_V_norm=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)
            local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}}
    T dt_curvature=(curvature_motion && T_GRID::dimension>1)?sigma*2*grid.one_over_dX.Magnitude_Squared():0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR FAST_LEVELSET<T_GRID>::
CFL(const T_ARRAYS_VECTOR& V) const
{
    T dt_convection=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi(cell))<=half_band_width)
        dt_convection=max(dt_convection,TV::Dot_Product(V(cell),grid.one_over_dX));}
    T dt_curvature=(curvature_motion && T_GRID::dimension>1)?sigma*2*grid.one_over_dX.Magnitude_Squared():0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Functions Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_LEVELSET<T_GRID>::
Fast_Marching_Method(const int local_advection_spatial_order, const T time)
{
    T_LEVELSET::Fast_Marching_Method(time,half_band_width+grid.dX.Max()*(1+min(3,local_advection_spatial_order)));
    boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented - pseudo-time
}
//#####################################################################
// Functions Fast_Marching_Method_Outside_Band
//#####################################################################
template<class T_GRID> void FAST_LEVELSET<T_GRID>::
Fast_Marching_Method_Outside_Band(const int local_advection_spatial_order, const T time)
{
    T_LEVELSET::Fast_Marching_Method_Outside_Band(half_band_width,time,half_band_width+grid.dX.Max()*(1+min(3,local_advection_spatial_order)));
    boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented - pseudo-time
}
//#####################################################################
template class FAST_LEVELSET<GRID<VECTOR<float,1> > >;
template class FAST_LEVELSET<GRID<VECTOR<float,2> > >;
template class FAST_LEVELSET<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_LEVELSET<GRID<VECTOR<double,1> > >;
template class FAST_LEVELSET<GRID<VECTOR<double,2> > >;
template class FAST_LEVELSET<GRID<VECTOR<double,3> > >;
#endif
