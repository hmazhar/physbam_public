//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_DYADIC
//#####################################################################
//
// Solves u_xx + u_yy = f.
// Input u as ARRAYS with values equal to the initial guess at Dirichlet boundary conditions points.
// In the Neumann case f must sum to 0 for compatibility.
// Set psi_D=true for each point that gets Dirichlet boundary conditions.
// Set psi_N=true for each face that gets Neumann boundary conditions.
//
//#####################################################################  
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __LAPLACE_DYADIC__
#define __LAPLACE_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T> class SPARSE_MATRIX_NXN;

template<class T_GRID>
class LAPLACE_DYADIC:public LAPLACE<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::INDEX INDEX;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
public:
    typedef LAPLACE<T> BASE;
    using BASE::tolerance;using BASE::number_of_regions;using BASE::solve_neumann_regions;

    const T_GRID& grid;
    ARRAY<T>& u;
    ARRAY<T> f;
    PCG_SPARSE<T> pcg;
    ARRAY<int> filled_region_colors;
    ARRAY<bool> filled_region_touches_dirichlet;
    ARRAY<bool> psi_N,psi_D;
    ARRAY<T> u_interface; // interface boundary condition - 2nd order method
    bool enforce_compatibility;

public:
    LAPLACE_DYADIC(const T_GRID& grid_input,ARRAY<T>& u_input,const bool initialize_grid,const bool enforce_compatibility_input)
        :grid(grid_input),u(u_input),enforce_compatibility(enforce_compatibility_input)
    {if(initialize_grid) Initialize_Grid();}

    virtual ~LAPLACE_DYADIC()
    {}

    void Initialize_Grid()
    {f.Resize(grid.number_of_cells);
    psi_N.Resize(grid.number_of_faces);psi_D.Resize(grid.number_of_cells);
    filled_region_colors.Resize(grid.number_of_cells);}
    
    bool All_Cell_Faces_Neumann(const CELL* cell)
    {for(int i=0;i<T_GRID::number_of_faces_per_cell;i++){
        ARRAY<CELL*> face_neighbors;cell->Get_All_Face_Neighbors(i,face_neighbors,&grid);
        for(int neighbor=1;neighbor<=face_neighbors.m;neighbor++){CELL* neighbor_cell=face_neighbors(neighbor);
            if(cell->Depth_Of_This_Cell()<=neighbor_cell->Depth_Of_This_Cell()){if(!psi_N(neighbor_cell->Face(MAP_MESH::opposite_face[i]))) return false;}
            else if(!psi_N(cell->Face(i))) return false;}}
    return true;}

//#####################################################################
    virtual void Solve(const T time,const bool solution_regions_already_computed=false);
    void Solve_Subregion(ARRAY<int>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b);
    virtual void Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index);
    virtual void Find_A_Part_One(const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts,ARRAY<ARRAY<T> >& row_sum);
    virtual void Find_A_Part_Two(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<T> >& row_sum);
public:
    void Find_Solution_Regions();
    void Set_Neumann_Outer_Boundaries();
    void Set_Dirichlet_Outer_Boundaries();
//#####################################################################
};
}
#endif
#endif
