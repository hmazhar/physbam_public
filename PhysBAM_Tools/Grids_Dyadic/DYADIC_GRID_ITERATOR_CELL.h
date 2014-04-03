//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_GRID_ITERATOR_CELL
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __DYADIC_GRID_ITERATOR_CELL__
#define __DYADIC_GRID_ITERATOR_CELL__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class TV> class RANGE;
template<class T_GRID>
class DYADIC_GRID_ITERATOR_CELL
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL CELL;typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;
    typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;
public:
    const T_GRID& grid;
    UNIFORM_CELL_ITERATOR uniform_cell_iterator;
    ARRAY<CELL*>& cell_pointer_from_index;
    ARRAY<VECTOR<CELL*,T_GRID::number_of_neighbors_per_cell> >& neighbors;
    CELL* current_cell;

    DYADIC_GRID_ITERATOR_CELL(const T_GRID& grid_input,const int number_of_ghost_cells=0,const typename UNIFORM_GRID::REGION& region_type=UNIFORM_GRID::WHOLE_REGION,const int side=0)
        :grid(grid_input),uniform_cell_iterator(grid_input.uniform_grid,number_of_ghost_cells,region_type,side),cell_pointer_from_index(grid.Cell_Pointer_From_Index()),neighbors(grid.Neighbors())
    {
        current_cell=Get_First_Cell_In_Tree(grid.cells(uniform_cell_iterator.Cell_Index()));
    }

    DYADIC_GRID_ITERATOR_CELL(const T_GRID& grid_input,const RANGE<TV_INT>& region)
        :grid(grid_input),uniform_cell_iterator(grid_input.uniform_grid,region),cell_pointer_from_index(grid.Cell_Pointer_From_Index()),neighbors(grid.Neighbors())
    {
        current_cell=Get_First_Cell_In_Tree(grid.cells(uniform_cell_iterator.Cell_Index()));
    }

    CELL* Get_First_Cell_In_Tree(CELL* cell) const
    {while(cell->Has_Children())cell=cell->Child(0);return cell;}

    void Next()
    {if(!current_cell->Parent()){
        uniform_cell_iterator.Next();
        if(uniform_cell_iterator.Valid())current_cell=Get_First_Cell_In_Tree(grid.cells(uniform_cell_iterator.Cell_Index()));
        return;}
    if(current_cell->child_index==T_GRID::number_of_children_per_cell-1){
        current_cell=current_cell->Parent();
        Next();return;}
    else{current_cell=Get_First_Cell_In_Tree(current_cell->Parent()->Child(current_cell->child_index+1));}}

    bool Valid() const
    {return uniform_cell_iterator.Valid();}

    int Cell_Index() const // returns a number between 1 and grid.number_of_cells. The global cell value
    {return current_cell->Cell();}

    CELL* Cell_Pointer() const
    {return current_cell;}

    CELL* Cell_Neighbor(const int i) const
    {return neighbors(Cell_Index())(i);}

    TV Location() const
    {return current_cell->Center();}

    TV DX() const
    {return current_cell->DX();}

//#####################################################################
};
}
#endif
#endif
