//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_FACE_HORIZONTAL
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_FACE_HORIZONTAL__
#define __RLE_GRID_ITERATOR_FACE_HORIZONTAL__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
namespace PhysBAM{

template<class T_GRID,int axis>
class RLE_GRID_ITERATOR_FACE_HORIZONTAL
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;
    typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;

    CELL_ITERATOR cell1,cell2;

    RLE_GRID_ITERATOR_FACE_HORIZONTAL(const T_GRID& grid,const int number_of_ghost_cells,const T_BOX_HORIZONTAL_INT& sentinels=T_BOX_HORIZONTAL_INT::Zero_Box())
        :cell1(grid,number_of_ghost_cells,sentinels-Sentinels()),cell2(grid,number_of_ghost_cells,sentinels+Sentinels())
    {}

    RLE_GRID_ITERATOR_FACE_HORIZONTAL(const T_GRID& grid,const int number_of_ghost_cells,const bool include_boundary)
        :cell1(grid,number_of_ghost_cells,Boundary_Sentinels(include_boundary)-Sentinels()),cell2(grid,number_of_ghost_cells,Boundary_Sentinels(include_boundary)+Sentinels())
    {}

    RLE_GRID_ITERATOR_FACE_HORIZONTAL(const T_GRID& grid,const T_BOX_HORIZONTAL_INT& region)
        :cell1(grid,region-Direction()),cell2(grid,region)
    {}

    static int Axis()
    {return axis;}

    static int Horizontal_Axis()
    {return axis/2+1;}

    static typename T_GRID::VECTOR_HORIZONTAL_INT Direction()
    {return TV_HORIZONTAL_INT::Axis_Vector(Horizontal_Axis());}

    static T_BOX_HORIZONTAL_INT Boundary_Sentinels(const bool include_boundary)
    {return include_boundary?T_BOX_HORIZONTAL_INT::Zero_Box():T_BOX_HORIZONTAL_INT(Direction(),-Direction());}

    static T_BOX_HORIZONTAL_INT Sentinels()
    {return T_BOX_HORIZONTAL_INT(TV_HORIZONTAL_INT(),Direction());}

    operator bool() const
    {assert((bool)cell1==(bool)cell2);return cell1;}

    void operator++(int)
    {if(cell1.jmax()==cell2.jmax()){cell1++;cell2++;}
    else if(cell1.jmax()<cell2.jmax()) cell1++;else cell2++;}

    int Face() const
    {return max(cell1.Face(axis,1),cell2.Face(axis,0));}

    bool Both_Cells_Short() const
    {return cell1.Short() && cell2.Short();}

    bool Cell_Short(const int side) const
    {assert(0<=side&&side<=1);return (side?cell2:cell1).Short();}

    int Cell(const int side) const
    {assert(0<=side&&side<=1);return side?cell2.Cell():cell1.Cell();}

    bool Long() const
    {return Length()>1;}

    bool Short() const
    {return !Long();}

    bool First_In_Column() const
    {return cell1.First_In_Column() && cell2.First_In_Column();}

    int i() const
    {return cell2.i;}

    int ij() const
    {return cell2.ij;}

    int j() const
    {return max(cell1.j,cell2.j);}

    int jmax() const
    {return min(cell1.jmax(),cell2.jmax());}

    int Length() const
    {return jmax()-j();}

    TV X() const
    {typename T_GRID::VECTOR_INT index=cell2.I();index.y=j();TV X=cell1.grid.uniform_grid.Face(axis,index);if(Long()) X.y+=(T).5*cell1.grid.uniform_grid.dX.y*(Length()-1);return X;}

    TV Cell_X(const int side) const
    {assert(Both_Cells_Short());assert(0<=side&&side<=1);return side?cell2.X():cell1.X();}

    static RANGE<VECTOR<int,1> > Indices_In_Slices(const T_GRID& grid,const int i_start,const int i_end)
    {assert(axis==1);return RANGE<VECTOR<int,1> >(grid.columns(i_start)(1).faces[0],grid.columns(i_end)(1).faces[1]-1);}

    static RANGE<VECTOR<int,1> > Indices_In_Slice(const T_GRID& grid,const int i,const int ij_start,const int ij_end)
    {int face1=2*(axis-1),ij=ij_end+(axis==1),face2=face1+(axis==3);
    return RANGE<VECTOR<int,1> >(grid.columns(i,ij_start)(1).faces[face1],grid.columns(i,ij)(1).faces[face2]-1);}

    static void Compute_Indices(T_GRID& grid,int& number_of_faces)
    {int offset=2*(axis-1);
    for(RLE_GRID_ITERATOR_FACE_HORIZONTAL<T_GRID,axis> face(grid,grid.number_of_ghost_cells);face;face++){number_of_faces++;
        if(!face.cell1.dj && face.cell1.j>=face.cell2.j) ((typename T_GRID::RUN*)face.cell1.run)->faces[offset+1]=number_of_faces;
        if(!face.cell2.dj && face.cell2.j>=face.cell1.j) ((typename T_GRID::RUN*)face.cell2.run)->faces[offset+0]=number_of_faces;
        if(face.Long()) number_of_faces+=grid.long_run_faces_horizontal-1;}}

//#####################################################################
};
}
#endif
#endif
