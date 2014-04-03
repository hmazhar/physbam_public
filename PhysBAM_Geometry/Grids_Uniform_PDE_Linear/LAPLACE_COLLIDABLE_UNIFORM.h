//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_UNIFORM
//#####################################################################
#ifndef __LAPLACE_COLLIDABLE_UNIFORM__
#define __LAPLACE_COLLIDABLE_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Grids_PDE_Linear/LAPLACE_COLLIDABLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class T_GRID>
class LAPLACE_COLLIDABLE_UNIFORM:public LAPLACE_UNIFORM<T_GRID>,public LAPLACE_COLLIDABLE<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL T_FLOOD_FILL;typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
public:
    typedef LAPLACE_UNIFORM<T_GRID> BASE;
    typedef LAPLACE_COLLIDABLE<T_GRID> COLLIDABLE_BASE;
    typedef TV VECTOR_T;
    typedef T_GRID GRID_T;

    using BASE::grid;using BASE::psi_N;
    using LAPLACE<T>::tolerance;using LAPLACE<T>::number_of_regions;using LAPLACE<T>::solve_neumann_regions;
    using COLLIDABLE_BASE::second_order_cut_cell_method;using COLLIDABLE_BASE::second_order_cut_cell_threshold;
    using COLLIDABLE_BASE::levelset;using COLLIDABLE_BASE::u_interface;

    //T_LEVELSET* levelset; // used in second order accurate cut cell method
    //T_FACE_ARRAYS_SCALAR u_interface; // interface boundary condition - 2nd order method
protected:
    T_ARRAYS_SCALAR phi_default;
    T_LEVELSET* levelset_default;
public:

    LAPLACE_COLLIDABLE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input,THREAD_QUEUE* thread_queue=0);
    LAPLACE_COLLIDABLE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,T_LEVELSET& cell_centered_levelset,const bool initialize_grid,const bool multiphase_input,
        const bool enforce_compatibility_input,THREAD_QUEUE* thread_queue=0);
    virtual ~LAPLACE_COLLIDABLE_UNIFORM();

    void Use_Internal_Level_Set()
    {levelset=levelset_default;phi_default.Resize(grid.Domain_Indices(1));}

    void Use_External_Level_Set(T_LEVELSET& cell_centered_levelset)
    {levelset=&cell_centered_levelset;phi_default.Clean_Memory();}

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid_input);
    virtual void Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input=true);
    virtual void Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
private:
    void Apply_Second_Order_Cut_Cell_Method(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
//#####################################################################
};
}
#endif
