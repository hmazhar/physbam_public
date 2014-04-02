//#####################################################################
// Copyright 2005, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_DYADIC
//#####################################################################
//
// Solves div (beta grad u) = f where beta > 0 with [u] and [beta un] given. 
// Input a LEVELSET_OCTREE class with phi as (1,number_of_nodes).
// Note that []=positive-negative according to phi.
//
//##################################################################### 
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __POISSON_DYADIC__
#define __POISSON_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic_PDE_Linear/LAPLACE_DYADIC.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/POISSON.h>
namespace PhysBAM{

template<class T_GRID>
class POISSON_DYADIC:public POISSON<typename T_GRID::SCALAR>,public LAPLACE_DYADIC<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::INDEX INDEX;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;
public:
    typedef T_GRID GRID_T;
    using POISSON<T>::u_jumps;using POISSON<T>::beta_minus;using POISSON<T>::beta_plus;using POISSON<T>::GFM;using POISSON<T>::smear_beta;using POISSON<T>::number_of_interface_cells;
    using POISSON<T>::beta_un_jumps;using POISSON<T>::use_variable_beta;using POISSON<T>::beta_given_on_faces;using LAPLACE_DYADIC<T_GRID>::psi_D;
    using LAPLACE_DYADIC<T_GRID>::psi_N;using LAPLACE_DYADIC<T_GRID>::f;using LAPLACE_DYADIC<T_GRID>::grid;
    using LAPLACE_DYADIC<T_GRID>::Find_Solution_Regions;using LAPLACE_DYADIC<T_GRID>::filled_region_colors;

    ARRAY<T> beta_face;
    ARRAY<T> u_jump,beta_un_jump; // [u] and [beta un] on the grid
    ARRAY<T> variable_beta;

protected:
    T dt;
    bool dt_is_set;

public:
    POISSON_DYADIC(const T_GRID& grid_input,ARRAY<T>& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input)
        :LAPLACE_DYADIC<T_GRID>(grid_input,u_input,false,enforce_compatibility_input)
    {
        Initialize_Grid(true);
    }
    
    virtual ~POISSON_DYADIC()
    {}
    
    void Set_Dt(const T dt_in)
    {dt=dt_in;dt_is_set=true;}

    void Set_Jump()
    {u_jump.Resize(grid.number_of_cells);u_jumps=true;}

    void Set_Derivative_Jump()
    {beta_un_jump.Resize(grid.number_of_cells);beta_un_jumps=true;}

    void Set_Variable_beta(const bool beta_given_on_faces_input=false)
    {use_variable_beta=true;beta_given_on_faces=beta_given_on_faces_input;
    if(!beta_given_on_faces)variable_beta.Resize(grid.number_of_cells);}

    void Initialize_Grid(const bool skip_laplace_grid_initialization=false) // skip flag for constructor because base class constructor already inits
    {if(!skip_laplace_grid_initialization)LAPLACE_DYADIC<T_GRID>::Initialize_Grid();
    beta_face.Resize(grid.number_of_faces);
    if(u_jumps)u_jump.Resize(grid.number_of_cells);
    if(beta_un_jumps)beta_un_jump.Resize(grid.number_of_cells);
    if(use_variable_beta)variable_beta.Resize(grid.number_of_cells);}
    
//#####################################################################
    void Compute_beta_And_Add_Jumps_To_b(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Find_Constant_beta();
    void Find_Variable_beta();
    virtual void Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index);
    virtual void Find_A_Part_One(const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts,ARRAY<ARRAY<T> >& row_sum);
    virtual void Find_A_Part_Two(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<T> >& row_sum);
//#####################################################################
};
}
#endif
#endif
