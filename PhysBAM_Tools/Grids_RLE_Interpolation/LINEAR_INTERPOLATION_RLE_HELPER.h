//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_RLE_HELPER
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LINEAR_INTERPOLATION_RLE_HELPER__
#define __LINEAR_INTERPOLATION_RLE_HELPER__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
namespace PhysBAM{

template<class T_GRID>
class LINEAR_INTERPOLATION_RLE_HELPER
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename T_GRID::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;

private:
    LINEAR_INTERPOLATION_RLE_HELPER(); // disallow construction
public:

    template<class T2>
    static T2 Interpolate_Cells_Short(const T_GRID& grid,const ARRAY<T2>& u,const TV& X)
    {T_BLOCK block(grid,X);assert(block);return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Cell(block,u,X);}

    // leaves long cell values undefined
    static void Interpolate_From_Faces_To_Short_Cells(const T_GRID& grid,const ARRAY<T>& u_face,ARRAY<TV>& u_cell)
    {for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();
        for(int axis=1;axis<=T_GRID::dimension;axis++)u_cell(c)[axis]=(T).5*(u_face(cell.First_Face_Index(axis))+u_face(cell.Second_Face_Index(axis)));}}

    // leaves long face values and boundary short face values undefined
    static void Interpolate_From_Short_Cells_To_Short_Faces(const T_GRID& grid,const ARRAY<T>& u_cell,ARRAY<T>& u_face)
    {T_GRID::template Face_Loop<Interpolate_From_Short_Cells_To_Short_Faces_Helper>(grid,u_cell,u_face);}

private:
    struct Interpolate_From_Short_Cells_To_Short_Faces_Helper{template<class T_FACE> static void Apply(const T_GRID& grid,const ARRAY<T>& u_cell,ARRAY<T>& u_face)
    {for(T_FACE face(grid,grid.number_of_ghost_cells,false);face;face++){int f=face.Face(),c1=face.cell1.Cell(),c2=face.cell2.Cell();
        if(face.Both_Cells_Short()) u_face(f)=(T).5*(u_cell(c1)+u_cell(c2));
        else if(face.cell1.Short()) u_face(f)=u_cell(c1);
        else if(face.cell2.Short()) u_face(f)=u_cell(c2);}}};

//#####################################################################
};
}
#endif
#endif
