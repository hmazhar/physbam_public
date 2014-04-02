//#####################################################################
// Copyright 2003-2005, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_DYADIC  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
// #if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __LEVELSET_DYADIC__
#define __LEVELSET_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic/BLOCK_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class T> class MATRIX_3X3;

template<class T_GRID_input>
class LEVELSET_DYADIC:public LEVELSET<typename T_GRID_input::SCALAR,T_GRID_input>
{
public:
    typedef T_GRID_input T_GRID;
private:
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;
    typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;typedef typename UNIFORM_GRID::NODE_ITERATOR UNIFORM_NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<UNIFORM_GRID>::ARRAYS_SCALAR UNIFORM_ARRAYS;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
public:
    typedef LEVELSET<T,T_GRID> BASE;
    using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;using BASE::small_number;using BASE::boundary;using BASE::max_time_step;
    using BASE::interpolation;using BASE::secondary_interpolation;
    using BASE::curvature_motion;
    using BASE::sigma;using BASE::collision_body_list;using BASE::collidable_phi_replacement_value;using BASE::valid_mask_current;
    
    T_GRID& grid;
    ARRAY<T>& phi;
    ARRAY<TV>* normals;
    ARRAY<T>* curvature;
    T half_band_width; // currently unused

    LEVELSET_DYADIC(T_GRID& grid_input,ARRAY<T>& phi_input);
    virtual ~LEVELSET_DYADIC();

    void Initialize_Levelset_Grid_Values() PHYSBAM_OVERRIDE
    {LEVELSET<T,T_GRID>::Initialize_Levelset_Grid_Values(grid);delete normals;normals=0;delete curvature;curvature=0;}
  
    void Set_Band_Width(const T number_of_cells=20) 
    {half_band_width=number_of_cells*grid.Minimum_Edge_Length()/2;}

    T Phi(const TV& location,const ARRAY<T>* phi_nodes=0) const 
    {return interpolation->Clamped_To_Array_Cell(grid,phi,phi_nodes,location);}

    T Phi(const CELL* cell,const TV& location,const ARRAY<T>& phi_nodes) const 
    {return interpolation->From_Leaf_Cell_Cell(grid,cell,phi,phi_nodes,location);}

    T Phi(const BLOCK_DYADIC<T_GRID>& block,const TV& location) const 
    {assert(block.Inside(location));return interpolation->From_Block_Cell(grid,block,phi,location);}

    T Phi_From_Close_Cell(const CELL* cell,const TV& location,const ARRAY<T>* phi_nodes=0) const 
    {return interpolation->From_Close_Cell_Cell(grid,cell,phi,phi_nodes,location);}

    T Phi_Secondary(const TV& location,const ARRAY<T>* phi_nodes=0) const
    {return secondary_interpolation->Clamped_To_Array_Cell(grid,phi,phi_nodes,location);}

    T Clamped_Phi(const TV& location,const ARRAY<T>* phi_nodes=0) const 
    {return interpolation->Clamped_To_Array_Cell(grid,phi,phi_nodes,grid.Clamp(location));}

    T Extended_Phi(const TV& location,const ARRAY<T>* phi_nodes=0) const
    {TV clamped_location(grid.Clamp(location));
    return interpolation->Clamped_To_Array_Cell(grid,phi,phi_nodes,clamped_location)+(clamped_location-location).Magnitude();}

    T Extended_Phi_From_Close_Cell(const CELL* cell,const TV& location,const ARRAY<T>* phi_nodes=0) const
    {TV clamped_location(grid.Clamp(location));
    return interpolation->From_Close_Cell_Cell(grid,cell,phi,phi_nodes,clamped_location)+(clamped_location-location).Magnitude();}

    bool Lazy_Inside(const TV& unclamped_X,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {return Phi(unclamped_X,phi_nodes)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {phi_value=Phi(unclamped_X,phi_nodes);return phi_value<=contour_value;}

    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {return Extended_Phi(unclamped_X,phi_nodes)<=contour_value;}

    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {phi_value=Extended_Phi(unclamped_X,phi_nodes);return phi_value<=contour_value;}

    bool Lazy_Outside(const TV& unclamped_X,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {return Phi(unclamped_X,phi_nodes)>contour_value;}

    bool Lazy_Outside_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {phi_value=Phi(unclamped_X,phi_nodes);return phi_value>contour_value;}

    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {return Extended_Phi(unclamped_X,phi_nodes)>contour_value;}

    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0,const ARRAY<T>* phi_nodes=0) const
    {phi_value=Extended_Phi(unclamped_X,phi_nodes);return phi_value>contour_value;}

//#####################################################################
    T Collision_Aware_Phi(const TV& location) const;
    T CFL(const ARRAY<T>& face_velocities) const;
    void Compute_Normals(const T time=0);
    TV Normal(const TV& location,const bool use_precomputed_normals_if_possible=true,const ARRAY<T>* phi_nodes=0) const;
    TV Normal_From_Leaf_Cell(const CELL* leaf_cell,const TV& location,const ARRAY<T>* phi_nodes=0) const;
    TV Normal(const BLOCK_DYADIC<T_GRID>& block,const TV& location,const bool use_precomputed_normals_if_possible=true,const ARRAY<T>* phi_nodes=0) const;
    TV Extended_Normal(const TV& location,const bool use_precomputed_normals_if_possible=true,const ARRAY<T>* phi_nodes=0) const;
    TV Extended_Normal_From_Leaf_Cell(const CELL* leaf_cell,const TV& location,const ARRAY<T>* phi_nodes=0) const;
    void Create_Dyadic_Levelset(const int max_depth,const int min_depth,T (*phi_sample_function)(const TV& location),const bool verbose=true);
    bool Whitney_Criteria(const CELL* hexahedron,const T half_band_width=0) const;
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<int>* seed_indices=0);
    void Get_Signed_Distance_Using_FMM(ARRAY<T>& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<int>* seed_indices=0);
    void Fast_Marching_Method_Outside_Band(const T half_band_width,const T time=0,const T stopping_distance=0);
public:
    T Approximate_Negative_Material(const T interface_thickness=3,const T time=0) const;
    T Approximate_Positive_Material(const T interface_thickness=3,const T time=0) const;
//#####################################################################
};
}
#endif
#endif
