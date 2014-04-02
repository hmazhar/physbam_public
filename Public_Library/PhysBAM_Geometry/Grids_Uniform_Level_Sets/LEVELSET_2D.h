//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_2D  
//##################################################################### 
#ifndef __LEVELSET_2D__
#define __LEVELSET_2D__ 

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID_input>
class LEVELSET_2D:public LEVELSET_UNIFORM<T_GRID_input>
{
    typedef typename T_GRID_input::SCALAR T;typedef typename T_GRID_input::CELL_ITERATOR CELL_ITERATOR;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef T_GRID_input T_GRID;    
    typedef LEVELSET_UNIFORM<T_GRID> BASE;
    using BASE::grid;using BASE::phi;using BASE::normals;using BASE::curvature;using BASE::cell_range;
    using BASE::collision_body_list;using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;
    using BASE::levelset_callbacks;using BASE::small_number;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature_motion;using BASE::sigma;using BASE::interpolation;
    using BASE::curvature_interpolation;using BASE::normal_interpolation;using BASE::collision_aware_interpolation_minus;using BASE::number_of_ghost_cells;

    LEVELSET_2D(T_GRID& grid_input,ARRAY<T,VECTOR<int,2> >& phi_input,const int number_of_ghost_cells_input=3);
    ~LEVELSET_2D();

    VECTOR<T,2> Normal(const VECTOR<T,2>& location) const
    {if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,location).Normalized();
    else return VECTOR<T,2>((Phi(VECTOR<T,2>(location.x+grid.dX.x,location.y))-Phi(VECTOR<T,2>(location.x-grid.dX.x,location.y)))/(2*grid.dX.x),
        (Phi(VECTOR<T,2>(location.x,location.y+grid.dX.y))-Phi(VECTOR<T,2>(location.x,location.y-grid.dX.y)))/(2*grid.dX.y)).Normalized();}

    VECTOR<T,2> Extended_Normal(const VECTOR<T,2>& location) const
    {if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,grid.Clamp(location)).Normalized();
    else return VECTOR<T,2>((Extended_Phi(VECTOR<T,2>(location.x+grid.dX.x,location.y))-Extended_Phi(VECTOR<T,2>(location.x-grid.dX.x,location.y)))/(2*grid.dX.x),
        (Extended_Phi(VECTOR<T,2>(location.x,location.y+grid.dX.y))-Extended_Phi(VECTOR<T,2>(location.x,location.y-grid.dX.y)))/(2*grid.dX.y)).Normalized();}

    static VECTOR<T,2> Normal_At_Node(const T_GRID& grid,const ARRAY<T,VECTOR<int,2> >& phi,const VECTOR<int,2>& index)
    {int i=index.x,j=index.y;return VECTOR<T,2>((phi(i+1,j)-phi(i-1,j))*grid.one_over_dX.x,(phi(i,j+1)-phi(i,j-1))*grid.one_over_dX.y).Normalized();}

//#####################################################################
    SYMMETRIC_MATRIX<T,2> Hessian(const VECTOR<T,2>& X) const;
    VECTOR<T,1> Principal_Curvatures(const VECTOR<T,2>& X) const;
    void Compute_Normals(const T time=0);
    void Compute_Curvature(const T time=0);
    T Compute_Curvature(const ARRAY<T,VECTOR<int,2> >& phi_input,const VECTOR<int,2>& index) const;
    T Compute_Curvature(const VECTOR<T,2>& location) const;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true);
public:
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,2> >* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,2> >& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,2> >* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Fast_Marching_Method_Outside_Band(const T half_band_width,const T time=0,const T stopping_distance=0);
public:
    T Approximate_Length(const T interface_thickness=3,const T time=0) const;
//#####################################################################
};
}
#endif
