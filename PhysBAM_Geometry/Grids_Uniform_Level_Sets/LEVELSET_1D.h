//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_1D  
//##################################################################### 
#ifndef __LEVELSET_1D__
#define __LEVELSET_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class LEVELSET_1D:public LEVELSET_UNIFORM<GRID<VECTOR<T_input,1> > >
{
    typedef T_input T;
    STATIC_ASSERT(IS_FLOATING_POINT<T>::value);
public:
    typedef VECTOR<T,1> TV;
    typedef LEVELSET_UNIFORM<GRID<TV> > BASE;
    using BASE::grid;using BASE::phi;using BASE::normals;using BASE::curvature;using BASE::cell_range;
    using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;
    using BASE::levelset_callbacks;using BASE::small_number;using BASE::boundary;
    using BASE::max_time_step;using BASE::number_of_ghost_cells;
    using BASE::interpolation;using BASE::curvature_interpolation;using BASE::normal_interpolation;

    LEVELSET_1D(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,1> >& phi_input,const int number_of_ghost_cells_input=3);
    ~LEVELSET_1D();

    VECTOR<T,1> Normal(const VECTOR<T,1>& location) const
    {if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,location).Normalized();
    else return VECTOR<T,1>((Phi(location+grid.dX)-Phi(location-grid.dX))/(2*grid.dX.x)).Normalized();}
        
    VECTOR<T,1> Extended_Normal(const VECTOR<T,1>& location) const
    {if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,grid.Clamp(location)).Normalized();
    else return VECTOR<T,1>((Extended_Phi(location+grid.dX)-Extended_Phi(location-grid.dX))/(2*grid.dX.x)).Normalized();}

    static VECTOR<T,1> Normal_At_Node(const GRID<TV>& grid,const ARRAY<T,VECTOR<int,1> >& phi,const VECTOR<int,1>& index)
    {int i=index.x;return VECTOR<T,1>((phi(i+1)-phi(i-1))*grid.one_over_dX.x).Normalized();}

    T Compute_Curvature(const ARRAY<T,VECTOR<int,1> >& phi_input,const VECTOR<int,1>& index) const
    {return 0;}

    T Compute_Curvature(const VECTOR<T,1>& location) const
    {return 0;}

//#####################################################################
    MATRIX<T,1> Hessian(const VECTOR<T,1>& X) const;
    VECTOR<T,0> Principal_Curvatures(const VECTOR<T,1>& X) const;
    void Compute_Normals(const T time=0);
    void Compute_Curvature(const T time=0);
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true);
public:
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,1> >* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,1> >& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,1> >* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false);
    void Fast_Marching_Method_Outside_Band(const T half_band_width,const T time=0,const T stopping_distance=0);
//#####################################################################
};
}
#endif
