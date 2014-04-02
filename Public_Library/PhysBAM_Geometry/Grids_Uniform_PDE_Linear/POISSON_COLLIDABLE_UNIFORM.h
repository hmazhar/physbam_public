//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_COLLIDABLE_UNIFORM
//#####################################################################
#ifndef __POISSON_COLLIDABLE_UNIFORM__
#define __POISSON_COLLIDABLE_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/POISSON_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_PDE_Linear/LAPLACE_COLLIDABLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class T_GRID>
//class POISSON_COLLIDABLE_UNIFORM:public POISSON<typename T_GRID::SCALAR>,public LAPLACE_COLLIDABLE_UNIFORM<T_GRID>
class POISSON_COLLIDABLE_UNIFORM:public POISSON_UNIFORM<T_GRID>,public LAPLACE_COLLIDABLE<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;

public:
    typedef POISSON_UNIFORM<T_GRID> BASE;
    typedef LAPLACE_COLLIDABLE<T_GRID> COLLIDABLE_BASE;
    //typedef LAPLACE_COLLIDABLE_UNIFORM<T_GRID> BASE;
    using POISSON<T>::GFM;using POISSON<T>::number_of_interface_cells;using POISSON<T>::smear_beta;using POISSON<T>::beta_minus;using POISSON<T>::beta_plus;using POISSON<T>::beta_multiphase;
    using POISSON<T>::u_jumps;using POISSON<T>::beta_un_jumps;using POISSON<T>::use_variable_beta;using POISSON<T>::beta_given_on_faces;
    using POISSON<T>::multiphase;

    using BASE::grid;using BASE::psi_N;using BASE::filled_region_colors;using BASE::beta_face;
    using BASE::Find_Variable_beta;
    
    using COLLIDABLE_BASE::second_order_cut_cell_method;using COLLIDABLE_BASE::second_order_cut_cell_threshold;
    using COLLIDABLE_BASE::levelset;using COLLIDABLE_BASE::u_interface;

    T_ARRAYS_SCALAR u_jump,beta_un_jump; // [u] and [beta un] on the grid
    T_FACE_ARRAYS_SCALAR beta_interface_face; // 2nd order method
    T_FACE_ARRAYS_SCALAR u_jump_face;
    //T_LEVELSET* levelset; // used in second order accurate cut cell method
    LEVELSET_MULTIPLE_UNIFORM<T_GRID>* levelset_multiple;
    //T_FACE_ARRAYS_SCALAR u_interface; // interface boundary condition - 2nd order method
private:
    ARRAY<T_ARRAYS_SCALAR> phis_default;
    LEVELSET_MULTIPLE_UNIFORM<T_GRID> levelset_multiple_default;
protected:
    T dt;
    bool dt_is_set;
public:

    POISSON_COLLIDABLE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input);
    POISSON_COLLIDABLE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,T_LEVELSET& cell_centered_levelset,const bool initialize_grid,const bool multiphase_input,
        const bool enforce_compatibility_input);
    virtual ~POISSON_COLLIDABLE_UNIFORM();

    void Set_Dt(const T dt_in)
    {dt=dt_in;dt_is_set=true;}

    void Set_Jump()
    {u_jump.Resize(grid.Domain_Indices(1),false,false);u_jump.Fill(0);u_jumps=true;}

    void Set_Jump_Multiphase()
    {u_jump_face.Resize(grid,1,false,false);u_jump_face.Fill(0);u_jumps=true;}

    void Set_Derivative_Jump()
    {beta_un_jump.Resize(grid.Domain_Indices(1));beta_un_jumps=true;}

    void Initialize_Grid()
    {} // TODO(jontg): Stub

//#####################################################################
    void Use_Internal_Level_Set(const int number_of_regions);
    void Update_Internal_Level_Set(LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple_input);
    void Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input=true) PHYSBAM_OVERRIDE;
    void Initialize_Grid(const T_GRID& grid_input) PHYSBAM_OVERRIDE;
    void Compute_beta_And_Add_Jumps_To_b(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Find_Constant_beta(T_FACE_ARRAYS_SCALAR& beta_face,const T_ARRAYS_SCALAR& phi_ghost);
    void Find_Constant_beta(const T_ARRAYS_SCALAR& phi_ghost);
    void Find_Constant_beta_Multiphase(ARRAY<T_ARRAYS_SCALAR>& phis_ghost);
    virtual void Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index) PHYSBAM_OVERRIDE;
private:
    void Add_Jump_To_b(const T_ARRAYS_SCALAR& phi_ghost);
    void Add_Jump_To_b_Multiphase(ARRAY<T_ARRAYS_SCALAR>& phis_ghost);
    void Add_Derivative_Jump_To_b(const T_ARRAYS_SCALAR& phi_ghost);
    void Apply_Second_Order_Cut_Cell_Method(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
//#####################################################################
};
}
#endif

