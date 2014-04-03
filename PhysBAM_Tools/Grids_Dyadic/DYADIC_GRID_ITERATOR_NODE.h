//#####################################################################
// Copyright 2005, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_GRID_ITERATOR_NODE 
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __DYADIC_GRID_ITERATOR_NODE__
#define __DYADIC_GRID_ITERATOR_NODE__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
namespace PhysBAM{

template<class T_GRID>
class DYADIC_GRID_ITERATOR_NODE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;
public:
    const T_GRID& grid;
    ARRAY<CELL*>& cell_pointer_from_index;
    ARRAY<VECTOR<CELL*,T_GRID::number_of_neighbors_per_cell> >& neighbors;
    ARRAY<int>* indirection_array;
    int indirection_index;
    int current_index;
    ARRAY<ARRAY<int>*> indirection_arrays;
    int current_indirection_array;

    DYADIC_GRID_ITERATOR_NODE(const T_GRID& grid_input,const ARRAY<ARRAY<int>*>& indirection_arrays_input)
        :grid(grid_input),cell_pointer_from_index(grid.Cell_Pointer_From_Index()),neighbors(grid.Neighbors())
    {
        grid.Node_Iterator_Data();    
        indirection_arrays=indirection_arrays_input;
        indirection_index=0;current_index=0;
        indirection_array=indirection_arrays(current_indirection_array=1);
        Next();
    }

    void Next()
    {if(indirection_index<indirection_array->m){current_index=(*indirection_array)(++indirection_index);return;}
    else if(++current_indirection_array>indirection_arrays.m)return;
    indirection_array=indirection_arrays(current_indirection_array);current_index=(*indirection_array)(indirection_index=1);}

    bool Valid() const
    {return current_indirection_array<=indirection_arrays.m;}

    CELL* Deepest_Cell() const
    {return cell_pointer_from_index(grid.node_iterator_deepest_cells(current_index));}

    int Deepest_Cell_Index() const
    {return grid.node_iterator_deepest_cells(current_index);}

    int Node_Index_In_Cell() const // returns a number between 0 and number_of_nodes_per_cell-1. Which of the nodes in the deepest cell is being iterated over
    {return grid.nodes(current_index);}

    int Node_Index() const // returns a number between 1 and grid.number_of_nodes. The global node value
    {return current_index;}

    TV Location() const
    {return grid.Node_Location(Node_Index_In_Cell(),Deepest_Cell());}

//#####################################################################
};
}
#endif
#endif
