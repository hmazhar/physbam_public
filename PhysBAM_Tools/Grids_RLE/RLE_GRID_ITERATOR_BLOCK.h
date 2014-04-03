//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_BLOCK
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_BLOCK__
#define __RLE_GRID_ITERATOR_BLOCK__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <cassert>
namespace PhysBAM{

template<class T_GRID>
class RLE_GRID_ITERATOR_BLOCK
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::BLOCK_RUN BLOCK_RUN;
    typedef typename T_GRID::BOX_HORIZONTAL_INT BOX_HORIZONTAL_INT;
public:
    const T_GRID& grid;
    const BLOCK_RUN *run,*run_end;
    int j_end,j,dj;

    RLE_GRID_ITERATOR_BLOCK(const T_GRID& grid_input,const TV_INT& I) // creates an iterator that cannot be incremented
        :grid(grid_input)
    {
        run=grid.Clamped_Block(I);if(!run){j_end=0;j=1;dj=0;return;}
        j_end=j=I.y;dj=j-run->jmin;
    }

protected:
    RLE_GRID_ITERATOR_BLOCK(const T_GRID& grid_input)
        :grid(grid_input)
    {}

    void Initialize(const int index,const int index_end)
    {run=&grid.block_runs(index);run_end=&grid.block_runs(index_end);
    j_end=run->jmax;j=run->jmin;dj=0;}

    void Initialize(const RLE_GRID_ITERATOR_BLOCK<T_GRID>& block) // creates an iterator that cannot be incremented
    {run=block.run;j_end=block.j_end;j=block.j;dj=block.dj;}
public:

    operator bool() const
    {return j<=j_end;}

    int Block() const
    {return run->block+dj;}

    int Cell(const int index) const
    {assert(*this && 0<=index&&index<T_GRID::number_of_cells_per_block);
    static const int lookup[8][2]={{0,0},{1,0},{0,1},{1,1},{2,0},{3,0},{2,1},{3,1}};
    return run->cells[lookup[index][0]]+lookup[index][1]+dj;}

    int Face_Y(const int index) const
    {assert(0<=index&&index<T_GRID::number_of_faces_per_block/T_GRID::dimension);
    static const int lookup[12][2]={{0,0},{1,0},{0,1},{1,1},{0,2},{1,2},{2,0},{3,0},{2,1},{3,1},{2,2},{3,2}};
    return run->cells[lookup[index][0]]+lookup[index][1]+dj;}

    int Face_X(const int index) const
    {assert(0<=index&&index<T_GRID::number_of_faces_per_block/T_GRID::dimension);
    static const int lookup[12][2]={{0,0},{1,0},{2,0},{0,1},{1,1},{2,1},{3,0},{4,0},{5,0},{3,1},{4,1},{5,1}};
    return run->faces_x[lookup[index][0]]+lookup[index][1]+dj;}

    int Face_Z(const int index) const
    {assert(0<=index&&index<T_GRID::number_of_faces_per_block/T_GRID::dimension&&T_GRID::dimension==3);
    static const int lookup[12][2]={{0,0},{1,0},{0,1},{1,1},{2,0},{3,0},{2,1},{3,1},{4,0},{5,0},{4,1},{5,1}};
    return run->faces_z[lookup[index][0]]+lookup[index][1]+dj;}

    // Functions for LINEAR_INTERPOLATION_MAC_HELPER
    TV One_Over_DX() const
    {return grid.uniform_grid.one_over_dX;}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_X_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {return face_value(1,Face_X(block_face_index));}

    template<class T2>
    T2 Face_X_Value(const ARRAY<T2>& face_value,const int block_face_index) const
    {return face_value(Face_X(block_face_index));}

    T& Face_X_Reference(ARRAY<T>& face_value,const int block_face_index) const
    {return face_value(1,Face_X(block_face_index));}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_Y_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {return face_value(2,Face_Y(block_face_index));}

    template<class T2>
    T2 Face_Y_Value(const ARRAY<T2>& face_value,const int block_face_index) const
    {return face_value(Face_Y(block_face_index));}

    T& Face_Y_Reference(ARRAY<T>& face_value,const int block_face_index) const
    {return face_value(2,Face_Y(block_face_index));}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_Z_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {assert(T_GRID::dimension==3);return face_value(3,Face_Z(block_face_index));}

    template<class T2>
    T2 Face_Z_Value(const ARRAY<T2>& face_value,const int block_face_index) const
    {assert(T_GRID::dimension==3);return face_value(Face_Z(block_face_index));}

    T& Face_Z_Reference(ARRAY<T>& face_value,const int block_face_index) const
    {assert(T_GRID::dimension==3);return face_value(3,Face_Z(block_face_index));}

//#####################################################################
};
}
#endif
#endif
