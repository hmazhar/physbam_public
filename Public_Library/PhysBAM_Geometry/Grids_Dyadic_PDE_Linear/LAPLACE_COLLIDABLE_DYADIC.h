//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_DYADIC
//#####################################################################
//
// Solves u_xx + u_yy = f.
// Input u as ARRAYS with values equal to the initial guess at Dirichlet boundary conditions points.
// In the Neumann case f must sum to 0 for compatibility.
// Set psi_D=true for each point that gets Dirichlet boundary conditions.
// Set psi_N=true for each face that gets Neumann boundary conditions.
//
//#####################################################################  
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LAPLACE_COLLIDABLE_DYADIC__
#define __LAPLACE_COLLIDABLE_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_NXN;

template<class T_GRID>
class LAPLACE_COLLIDABLE_DYADIC:public LAPLACE<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
public:
    typedef LAPLACE<T> BASE;
    using BASE::tolerance;using BASE::number_of_regions;using BASE::solve_neumann_regions;

    T_GRID& grid;
    ARRAY<T>& u;
    ARRAY<T> f;
    PCG_SPARSE<T> pcg;
    ARRAY<int> filled_region_colors;
    ARRAY<bool> filled_region_touches_dirichlet;
    ARRAY<bool> psi_N,psi_D;
    T_LEVELSET* levelset;
    ARRAY<T> u_interface; // interface boundary condition - 2nd order method
    bool use_second_order_pressure;
protected:
    ARRAY<T> phi_default;
    T_LEVELSET levelset_default;

public:
    bool second_order_cut_cell_method;
    T second_order_cut_cell_threshold;

    LAPLACE_COLLIDABLE_DYADIC(T_GRID& grid_input,ARRAY<T>& u_input)
        :grid(grid_input),u(u_input),levelset(0),use_second_order_pressure(true),levelset_default(grid,phi_default),
        second_order_cut_cell_method(false),second_order_cut_cell_threshold((T)1e-3)
    {}

    virtual ~LAPLACE_COLLIDABLE_DYADIC()
    {}

    void Initialize_Grid()
    {f.Resize(grid.number_of_cells);
    psi_N.Resize(grid.number_of_faces);psi_D.Resize(grid.number_of_cells);
    filled_region_colors.Resize(grid.number_of_cells);
    if(levelset == &levelset_default)phi_default.Resize(grid.number_of_cells);
    Set_Up_Second_Order_Cut_Cell_Method(second_order_cut_cell_method);}
    
    void Use_Internal_Level_Set()
    {levelset=&levelset_default;phi_default.Resize(grid.number_of_cells);}

    void Use_External_Level_Set(T_LEVELSET& levelset_input)
    {levelset=&levelset_input;}

    void Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input=true)
    {second_order_cut_cell_method=use_second_order_cut_cell_method_input;
    if(second_order_cut_cell_method) u_interface.Resize(grid.number_of_faces);
    else u_interface.Resize(0);}

    bool All_Cell_Faces_Neumann(const CELL* cell)
    {for(int i=0;i<T_GRID::number_of_faces_per_cell;i++){
        ARRAY<CELL*> face_neighbors;cell->Get_All_Face_Neighbors(i,face_neighbors,&grid);
        for(int neighbor=1;neighbor<=face_neighbors.m;neighbor++){CELL* neighbor_cell=face_neighbors(neighbor);
            if(cell->Depth_Of_This_Cell()<=neighbor_cell->Depth_Of_This_Cell()){if(!psi_N(neighbor_cell->Face(MAP_MESH::opposite_face[i]))) return false;}
            else if(!psi_N(cell->Face(i))) return false;}}
    return true;}

//#####################################################################
    virtual void Solve(const T time,const bool solution_regions_already_computed=false);
    void Solve_Subregion(ARRAY<int>& matrix_index_to_cell_index,SPARSE_MATRIX_NXN<T>& A,VECTOR_ND<T>& b);
    virtual void Find_A(ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index);
private:
    void Find_A_Off_Diagonal_Helper(FACE_ITERATOR& iterator,ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,
        ARRAY<ARRAY<T> >& row_sum_array);
public:
    void Find_Solution_Regions();
    void Set_Neumann_Outer_Boundaries();
    void Set_Dirichlet_Outer_Boundaries();
//#####################################################################
};
}
#endif
#endif
