//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BLOCK_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __BLOCK_DYADIC__
#define __BLOCK_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <cassert>
namespace PhysBAM{

template<class TV> class RANGE;
template<class T_GRID>
class BLOCK_DYADIC
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL T_CELL;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    const T_GRID& grid;
    const T_CELL* cells[T_GRID::number_of_cells_per_block];

    BLOCK_DYADIC(const T_GRID& grid_input,const T_CELL* base_cell)
        :grid(grid_input)
    {
        Initialize(base_cell);
    }

    BLOCK_DYADIC(const T_GRID& grid_input,const TV& X)
        :grid(grid_input)
    {
        Initialize(grid.Base_Cell(X));
    }

    int Block() const
    {return cells[0]->Cell();}

    void Initialize(const T_CELL* base_cell)
    {assert(base_cell && grid.fully_refined_block(base_cell->Cell()));
    grid.All_Cells_In_Cell_Block(base_cell,cells);}

private:
    int Face_X(const int face_index) const
    {assert(0<=face_index&&face_index<T_GRID::number_of_faces_per_block/T_GRID::dimension);
    static const int lookup[][2]={{0,0},{0,1},{1,1},{2,0},{2,1},{3,1},{4,0},{4,1},{5,1},{6,0},{6,1},{7,1}};
    return cells[lookup[face_index][0]]->Face(lookup[face_index][1]);}

    int Face_Y(const int face_index) const
    {assert(0<=face_index&&face_index<T_GRID::number_of_faces_per_block/T_GRID::dimension);
    static const int lookup[][2]={{0,2},{1,2},{0,3},{1,3},{2,3},{3,3},{4,2},{5,2},{4,3},{5,3},{6,3},{7,3}};
    return cells[lookup[face_index][0]]->Face(lookup[face_index][1]);}

    int Face_Z(const int face_index) const
    {assert(0<=face_index&&face_index<T_GRID::number_of_faces_per_block/T_GRID::dimension&&T_GRID::dimension==3);
    static const int lookup[][2]={{0,4},{1,4},{2,4},{3,4},{0,5},{1,5},{2,5},{3,5},{4,5},{5,5},{6,5},{7,5}};
    return cells[lookup[face_index][0]]->Face(lookup[face_index][1]);}

public:
    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_X_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {return face_value(Face_X(block_face_index));}

    T& Face_X_Reference(T_FACE_ARRAYS_SCALAR& face_value,const int block_face_index) const
    {return face_value(Face_X(block_face_index));}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_Y_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {return face_value(Face_Y(block_face_index));}

    T& Face_Y_Reference(T_FACE_ARRAYS_SCALAR& face_value,const int block_face_index) const
    {return face_value(Face_Y(block_face_index));}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_Z_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {assert(T_GRID::dimension==3);return face_value(Face_Z(block_face_index));}

    T& Face_Z_Reference(T_FACE_ARRAYS_SCALAR& face_value,const int block_face_index) const
    {assert(T_GRID::dimension==3);return face_value(Face_Z(block_face_index));}

    int Cell(const int cell_index) const
    {assert(0<=cell_index&&cell_index<T_GRID::number_of_cells_per_block);
    return cells[cell_index]->Cell();}

    void All_Cell_Indices(int cell_indices[T_GRID::number_of_cells_per_block]) const
    {for(int i=0;i<T_GRID::number_of_cells_per_block;i++)cell_indices[i]=cells[i]->Cell();}

    const T_CELL* Base_Cell() const
    {return cells[0];}

    TV Minimum_Corner() const
    {return cells[0]->Center();}

    TV Center() const
    {return grid.Node_Location(T_GRID::number_of_nodes_per_cell-1,cells[0]);}

    TV One_Over_DX() const
    {return grid.One_Over_Minimum_Cell_DX();}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(Minimum_Corner(),cells[T_GRID::number_of_cells_per_block-1]->Center());}

    // TODO: these are not the typical parameters to Inside functions
    bool Inside(const TV& X,const T thickness_multiplier=-(T)1e-3) const
    {return Bounding_Box().Inside(X,thickness_multiplier*grid.Minimum_Edge_Length());}

    bool Lazy_Inside(const TV& X) const
    {return Bounding_Box().Lazy_Inside(X);}

//#####################################################################
};
}
#endif
#endif
