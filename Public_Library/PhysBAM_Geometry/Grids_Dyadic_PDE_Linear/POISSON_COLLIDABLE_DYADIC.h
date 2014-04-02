//#####################################################################
// Copyright 2005, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_COLLIDABLE_DYADIC
//#####################################################################
//
// Solves div (beta grad u) = f where beta > 0 with [u] and [beta un] given. 
// Input a LEVELSET_OCTREE class with phi as (1,number_of_nodes).
// Note that []=positive-negative according to phi.
//
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __POISSON_COLLIDABLE_DYADIC__
#define __POISSON_COLLIDABLE_DYADIC__

#include <PhysBAM_Tools/Grids_PDE_Linear/POISSON.h>
#include <PhysBAM_Geometry/Grids_Dyadic_PDE_Linear/LAPLACE_COLLIDABLE_DYADIC.h>
namespace PhysBAM{

template<class T_GRID>
class POISSON_COLLIDABLE_DYADIC:public POISSON<typename T_GRID::SCALAR>,public LAPLACE_COLLIDABLE_DYADIC<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
public:
    using POISSON<T>::u_jumps;using POISSON<T>::beta_minus;using POISSON<T>::beta_plus;using POISSON<T>::GFM;using POISSON<T>::smear_beta;using POISSON<T>::number_of_interface_cells;
    using POISSON<T>::beta_un_jumps;using POISSON<T>::use_variable_beta;using POISSON<T>::beta_given_on_faces;using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::psi_D;
    using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::psi_N;using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::f;using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::grid;
    using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::Find_Solution_Regions;using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::levelset;
    using LAPLACE_COLLIDABLE_DYADIC<T_GRID>::Use_Internal_Level_Set;

    ARRAY<T> beta_face;
    ARRAY<T> u_jump,beta_un_jump; // [u] and [beta un] on the grid
    ARRAY<T> variable_beta;

    POISSON_COLLIDABLE_DYADIC(T_GRID& grid_input,ARRAY<T>& u_input)
        :LAPLACE_COLLIDABLE_DYADIC<T_GRID>(grid_input,u_input)
    {
        Initialize_Grid(true);
    }
    
    virtual ~POISSON_COLLIDABLE_DYADIC()
    {}
    
    void Set_Jump()
    {u_jump.Resize(grid.number_of_cells);u_jumps=true;}

    void Set_Derivative_Jump()
    {beta_un_jump.Resize(grid.number_of_cells);beta_un_jumps=true;}

    void Set_Variable_beta(const bool beta_given_on_faces_input=false)
    {use_variable_beta=true;beta_given_on_faces=beta_given_on_faces_input;
    if(!beta_given_on_faces)variable_beta.Resize(grid.number_of_cells);}

    void Set_Up_Second_Order_Cut_Cell_Method()
    {PHYSBAM_NOT_IMPLEMENTED();} // not supported
    
    void Initialize_Grid(const bool skip_laplace_grid_initialization=false) // skip flag for constructor because base class constructor already inits
    {if(!skip_laplace_grid_initialization)LAPLACE_COLLIDABLE_DYADIC<T_GRID>::Initialize_Grid();
    beta_face.Resize(grid.number_of_faces);
    if(u_jumps)u_jump.Resize(grid.number_of_cells);
    if(beta_un_jumps)beta_un_jump.Resize(grid.number_of_cells);
    if(use_variable_beta)variable_beta.Resize(grid.number_of_cells);}
    
//#####################################################################
    void Compute_beta_And_Add_Jumps_To_b(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Find_Constant_beta(const ARRAY<T>& phi_ghost);
    void Find_Variable_beta();
    void Find_A(ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index) PHYSBAM_OVERRIDE;
private:
    void Find_A_Off_Diagonal_Helper(FACE_ITERATOR& iterator,ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,
        ARRAY<ARRAY<T> >& row_sum_array);
    void Add_Jump_To_b(const ARRAY<T>& phi_ghost);
    void Add_Derivative_Jump_To_b(const ARRAY<T>& phi_ghost);
//#####################################################################
};
}
#endif

#endif
