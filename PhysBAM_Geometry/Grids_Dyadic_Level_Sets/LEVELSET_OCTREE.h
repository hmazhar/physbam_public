//#####################################################################
// Copyright 2003-2005, Zhaosheng Bao, Ron Fedkiw, Sergey Koltakov, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_OCTREE  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LEVELSET_OCTREE__
#define __LEVELSET_OCTREE__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_DYADIC.h>
namespace PhysBAM{

template<class T,int d> class SYMMETRIC_MATRIX;

template<class T_input>
class LEVELSET_OCTREE:public LEVELSET_DYADIC<OCTREE_GRID<T_input> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef LEVELSET_DYADIC<OCTREE_GRID<T> > BASE;
protected:
    using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;
public:
    typedef OCTREE_GRID<T> T_GRID;
    using BASE::grid;using BASE::normals;using BASE::phi;using BASE::curvature;using BASE::interpolation;using BASE::boundary;using BASE::small_number;
    using BASE::collision_body_list;using BASE::clamp_phi_with_collision_bodies;

    LEVELSET_OCTREE(OCTREE_GRID<T>& grid_input,ARRAY<T>& phi_input)
        :LEVELSET_DYADIC<OCTREE_GRID<T> >(grid_input,phi_input)
    {}

    ~LEVELSET_OCTREE()
    {}

//#####################################################################
    SYMMETRIC_MATRIX<T,3> Hessian(const TV& X,const ARRAY<T>* phi_nodes=0) const;
    VECTOR<T,2> Principal_Curvatures(const TV& X,const ARRAY<T>* phi_nodes=0) const;
    void Compute_Curvature(const T time=0);
    T Curvature(const TV& location,const bool use_precomputed_curvatures_if_possible=true,const ARRAY<T>* phi_nodes=0) const;
//#####################################################################
};
}
#endif
#endif
