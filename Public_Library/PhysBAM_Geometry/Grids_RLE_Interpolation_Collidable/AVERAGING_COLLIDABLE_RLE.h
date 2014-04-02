//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __AVERAGING_COLLIDABLE_RLE__
#define __AVERAGING_COLLIDABLE_RLE__

#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP>
class AVERAGING_COLLIDABLE_RLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::BLOCK T_BLOCK;
public:
    typedef T_FACE_LOOKUP FACE_LOOKUP;

    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    T default_replacement_value;

    AVERAGING_COLLIDABLE_RLE(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const T default_replacement_value_input)
        :body_list(body_list_input),default_replacement_value(default_replacement_value_input)
    {}

    ~AVERAGING_COLLIDABLE_RLE()
    {}

    TV Face_To_Cell_Vector(const T_BLOCK& block,const CELL_ITERATOR& cell,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell);lookup.Set_Reference_Point(block,cell.X());
    TV value;for(int axis=1;axis<=T_GRID::dimension;axis++)
        value[axis]=(T).5*(lookup(axis,cell.First_Face_Index(axis)))+lookup(axis,cell.Second_Face_Index(axis));
    return value;}

/*
    T Cell_To_Face(const T_GRID& grid,const FACE_ITERATOR& face,const T_ARRAYS_SCALAR& u_cell) const // this never needs to set starting points doesn't use velocities
    {assert(iterator.Both_Cells_Short());
    if(body_list.Occupied_Face_Center(iterator)){
        T cell1_value=default_replacement_value,cell2_value=default_replacement_value;
        if(body_list.Cells_Center_Visible_From_Face(face.cell1,face.axis,face_index)) cell1_value=u_cell(face.cell1.Cell());
        if(body_list.Cells_Center_Visible_From_Face(face.cell2,face.axis,face_index)) cell2_value=u_cell(face.cell2.Cell());
        return (T).5*(cell1_value+cell2_value);}
    else return (T).5*(u_cell(face.cell1.Cell())+u_cell(face.cell2.Cell()));}
*/

    template<class T_FACE>
    TV Face_To_Face_Vector(const T_BLOCK& block,const T_FACE& face,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(face);lookup.Set_Reference_Point(block,face.X());
    return AVERAGING_RLE<T_GRID,T_FACE_LOOKUP>::Average_Face_To_Face_Vector_Helper(face,lookup,typename T_GRID::GRID_TAG());}

//#####################################################################
};
}
#endif
#endif
