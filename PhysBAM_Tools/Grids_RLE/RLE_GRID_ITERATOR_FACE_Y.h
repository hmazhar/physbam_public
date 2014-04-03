//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_FACE_Y
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_FACE_Y__
#define __RLE_GRID_ITERATOR_FACE_Y__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
namespace PhysBAM{

template<class T_GRID>
class RLE_GRID_ITERATOR_FACE_Y
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::RUN T_RUN;typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;
    typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;

    CELL_ITERATOR cell1,cell2;

    RLE_GRID_ITERATOR_FACE_Y(const T_GRID& grid,const int number_of_ghost_cells,const T_BOX_HORIZONTAL_INT& sentinels)
        :cell1(grid,number_of_ghost_cells,sentinels,-1,0),cell2(grid,number_of_ghost_cells,sentinels,0,1)
    {}

    RLE_GRID_ITERATOR_FACE_Y(const T_GRID& grid,const int number_of_ghost_cells,const bool include_boundary=false)
        :cell1(grid,number_of_ghost_cells,-Sentinels(),0-include_boundary,include_boundary-1),cell2(grid,number_of_ghost_cells,Sentinels(),1-include_boundary,include_boundary)
    {}

    RLE_GRID_ITERATOR_FACE_Y(const T_GRID& grid,const T_BOX_HORIZONTAL_INT& region,const bool include_boundary=false)
        :cell1(grid,region,0-include_boundary,include_boundary-1),cell2(grid,region,1-include_boundary,include_boundary)
    {}

    static T_BOX_HORIZONTAL_INT Sentinels()
    {return TV_HORIZONTAL_INT();}

    static int Axis()
    {return 2;}

    operator bool() const
    {assert((bool)cell1==(bool)cell2);return cell1;}

    void operator++(int)
    {cell1++;cell2++;assert(cell1.i==cell2.i);}

    int Face() const
    {return cell2.Face_Y();}

    int Cell(const int side=1) const
    {assert(0<=side&&side<=1);return cell2.Cell()+side-1;}

    bool Both_Cells_Short() const
    {return cell2.dj>0;}

    bool Lower_Cell_Long() const
    {return !cell2.First_In_Column() && cell1.Long();}

    bool Upper_Cell_Long() const
    {return !cell2.Last_In_Column() && cell2.Long();}

    bool Long() const
    {return Upper_Cell_Long();}

    int i() const
    {return cell2.i;}

    int ij() const
    {return cell2.ij;}

    int j() const
    {return cell2.j;}

    int jmax() const
    {return cell2.jmax();}

    int Length() const // for templatization purposes
    {return 1;}

    TV X() const
    {return cell2.grid.uniform_grid.Y_Face(cell2.I());}

    TV Cell_X(const int side) const
    {assert(Both_Cells_Short());assert(0<=side&&side<=1);return side?cell2.X():cell1.X();}

    static RANGE<VECTOR<int,1> > Indices_In_Slices(const T_GRID& grid,const int i_start,const int i_end)
    {return grid.Cells_In_Slices(i_start,i_end);}

    static RANGE<VECTOR<int,1> > Indices_In_Slice(const T_GRID& grid,const int i,const int ij_start,const int ij_end)
    {return grid.Cells_In_Slice(i,ij_start,ij_end);}

    static int Compute_Indices(T_GRID& grid)
    {int number_of_faces_y=0;
    for(RLE_GRID_ITERATOR_FACE_Y<T_GRID> face(grid,grid.number_of_ghost_cells,true);face;face++){number_of_faces_y++;
        if(!face.cell2.dj){
            T_RUN& run=*(T_RUN*)face.cell2.run;
            run.cell=number_of_faces_y; // face numbers are the same as cell numbers
            run.faces[2]=run.cell;run.faces[3]=run.cell+1;}
        if(face.Upper_Cell_Long()) number_of_faces_y+=grid.long_run_cells-1;}
    return number_of_faces_y;}

//#####################################################################
};
}
#endif
#endif
