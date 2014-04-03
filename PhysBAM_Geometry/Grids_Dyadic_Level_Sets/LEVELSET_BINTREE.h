//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_BINTREE  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
// #if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __LEVELSET_BINTREE__
#define __LEVELSET_BINTREE__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_DYADIC.h>
namespace PhysBAM{

template<class T,int d> class SYMMETRIC_MATRIX;

template<class T_input>
class LEVELSET_BINTREE:public LEVELSET_DYADIC<BINTREE_GRID<T_input> >
{
    typedef T_input T;
    typedef VECTOR<T,1> TV;
    typedef LEVELSET_DYADIC<BINTREE_GRID<T> > BASE;
protected:
    using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;
public:
    typedef BINTREE_GRID<T> T_GRID;
    using BASE::grid;using BASE::normals;using BASE::phi;using BASE::curvature;using BASE::interpolation;using BASE::boundary;using BASE::small_number;
    using BASE::collision_body_list;using BASE::clamp_phi_with_collision_bodies;

    LEVELSET_BINTREE(BINTREE_GRID<T>& grid_input,ARRAY<T>& phi_input)
        :LEVELSET_DYADIC<BINTREE_GRID<T> >(grid_input,phi_input)
    {}

    ~LEVELSET_BINTREE()
    {}

    SYMMETRIC_MATRIX<T,1> Hessian(const TV& X,const ARRAY<T>* phi_nodes=0) const
    {
        TV DX(grid.Leaf_Cell(X)->DX());
        T one_over_dx=1/DX.x;T two_phi_center=2*Phi(X,phi_nodes);
        T phi_xx=(Phi(TV(X.x+DX.x),phi_nodes)-two_phi_center+Phi(TV(X.x-DX.x),phi_nodes))*sqr(one_over_dx);
        return SYMMETRIC_MATRIX<T,1>(phi_xx);
    }

//#####################################################################
    VECTOR<T,0> Principal_Curvatures(const TV& X,const ARRAY<T>* phi_nodes=0) const {return VECTOR<T,0>();}
    void Compute_Curvature(const T time=0) {};
    T Curvature(const TV& location,const bool use_precomputed_curvatures_if_possible=true,const ARRAY<T>* phi_nodes=0) const {return 0;}
//#####################################################################
};
}
#endif
#endif
