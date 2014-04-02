//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_RLE
//#####################################################################
//
// Solves -A laplace u = -f, where A is the area of each cell
// Input u with values equal to the initial guess in Dirchlet boundary condition cells.
// In the Neumann case f must sum to 0 for compatibility.
// Set psi_D=true for each cell that gets Dirichlet boundary conditions.
// Set psi_N=true for each face that gets Neumann boundary conditions.
//
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LAPLACE_COLLIDABLE_RLE__
#define __LAPLACE_COLLIDABLE_RLE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Geometry/Parallel_Computation/LAPLACE_COLLIDABLE_RLE_MPI.h>
namespace PhysBAM{

class GRAPH;
template<class T_GRID> class LAPLACE_COLLIDABLE_RLE_MPI;
template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T_GRID> class MPI_RLE_GRID;

template<class T_GRID>
class LAPLACE_COLLIDABLE_RLE:public LAPLACE<typename T_GRID::SCALAR>
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;

    using LAPLACE<T>::tolerance;using LAPLACE<T>::number_of_regions;using LAPLACE<T>::solve_neumann_regions;

    const T_GRID& grid;
    ARRAY<T>& u;
    ARRAY<T> f;
    PCG_SPARSE<T> pcg;
    ARRAY<int> filled_region_colors;
    ARRAY<bool> filled_region_touches_dirichlet;
    ARRAY<bool> psi_N,psi_D;
    ARRAY<T> u_interface;
    LAPLACE_COLLIDABLE_RLE_MPI<T_GRID>* laplace_mpi;
    MPI_RLE_GRID<T_GRID>* mpi_grid;
    bool second_order_cut_cell_method;
    T second_order_cut_cell_threshold;

    LAPLACE_COLLIDABLE_RLE(const T_GRID& grid_input,ARRAY<T>& u_input)
        :grid(grid_input),u(u_input),mpi_grid(0),second_order_cut_cell_method(false),second_order_cut_cell_threshold((T)1e-3)
    {
        laplace_mpi=new LAPLACE_COLLIDABLE_RLE_MPI<T_GRID>(*this);
    }

    virtual ~LAPLACE_COLLIDABLE_RLE()
    {delete laplace_mpi;}

    virtual void Initialize_Grid()
    {f.Resize(grid.number_of_cells);psi_D.Resize(grid.number_of_cells,false,false);psi_N.Resize(grid.number_of_faces,false,false);
    filled_region_colors.Resize(grid.number_of_cells,false,false);Set_Up_Second_Order_Cut_Cell_Method(second_order_cut_cell_method);}

    void Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input=true)
    {second_order_cut_cell_method=use_second_order_cut_cell_method_input;
    u_interface.Resize(second_order_cut_cell_method?grid.number_of_faces:0);}

//#####################################################################
    void Solve(const T time=0,const bool solution_regions_already_computed=false,const ARRAY<T>* phi_ghost=0);
    void Find_Solution_Regions();
private:
    void Solve_Subregion(const ARRAY<int>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int color);
    void Find_Row_Lengths(ARRAY<ARRAY<int> >& row_lengths_array,ARRAY<int>& cell_index_to_matrix_index);
    struct Horizontal_Row_Length_Contribution{template<class T_FACE> static void Apply(const LAPLACE_COLLIDABLE_RLE<T_GRID>& laplace,ARRAY<ARRAY<int> >& row_lengths_array,
        ARRAY<int>& cell_index_to_matrix_index);};
    void Vertical_Row_Length_Contribution(ARRAY<ARRAY<int> >& row_lengths_array,ARRAY<int>& cell_index_to_matrix_index) const;
    void Find_A(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,const ARRAY<T>* phi_ghost);
    struct Find_Horizontal_Terms{template<class T_FACE> static void Apply(const LAPLACE_COLLIDABLE_RLE<T_GRID>& laplace,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,
        ARRAY<int>& cell_index_to_matrix_index,const ARRAY<T>* phi_ghost);};
    void Find_Vertical_Terms(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,const ARRAY<T>* phi_ghost) const;
    struct Add_Horizontal_Edges{template<class T_FACE> static void Apply(const LAPLACE_COLLIDABLE_RLE<T_GRID>& laplace,GRAPH& graph);};
//#####################################################################
};
}
#endif
#endif
