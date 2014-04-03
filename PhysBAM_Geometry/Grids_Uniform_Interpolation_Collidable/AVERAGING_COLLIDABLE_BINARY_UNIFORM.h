//#####################################################################
// Copyright 2009, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __AVERAGING_COLLIDABLE_BINARY_UNIFORM__
#define __AVERAGING_COLLIDABLE_BINARY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP>
class AVERAGING_COLLIDABLE_BINARY_UNIFORM
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:
    typedef T_FACE_LOOKUP FACE_LOOKUP;

    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    T default_replacement_value;

    AVERAGING_COLLIDABLE_BINARY_UNIFORM(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const T default_replacement_value_input)
        :body_list(body_list_input),default_replacement_value(default_replacement_value_input)
    {}

    ~AVERAGING_COLLIDABLE_BINARY_UNIFORM()
    {}

    TV Face_To_Node_Vector(const T_GRID& grid,const TV_INT& node_index,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(node_index);
    TV value;for(int axis=1;axis<=T_GRID::dimension;axis++)for(int face=1;face<=T_GRID::number_of_nodes_per_face;face++)
        value[axis]+=lookup(axis,T_GRID::Node_Face_Index(axis,node_index,face));
    return value/T_GRID::number_of_nodes_per_face;}

    TV Face_To_Cell_Vector(const T_GRID& grid,const TV_INT& cell_index,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell_index);lookup.Set_Reference_Point(grid.X(cell_index));
    TV value;for(int axis=1;axis<=T_GRID::dimension;axis++)
        value[axis]=(T).5*(lookup(axis,grid.First_Face_Index_In_Cell(axis,cell_index))+lookup(axis,grid.Second_Face_Index_In_Cell(axis,cell_index)));
    return value;}

    T Cell_To_Face(const T_GRID& grid,const int axis,const TV_INT& face_index,const T_ARRAYS_SCALAR& u_cell) const // this never needs to set starting points doesn't use velocities
    {FACE_ITERATOR iterator(grid,axis,face_index);
    TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face_index,cell1,cell2);
    if(body_list.Occupied_Face_Center(iterator)){
        T cell1_value=default_replacement_value,cell2_value=default_replacement_value;
        if(body_list.Cell_Center_Visible_From_Face(cell1,axis,face_index)) cell1_value=u_cell(cell1);
        if(body_list.Cell_Center_Visible_From_Face(cell2,axis,face_index)) cell2_value=u_cell(cell2);
        return (T).5*(cell1_value+cell2_value);}
    else return (T).5*(u_cell(cell1)+u_cell(cell2));}

    T Cell_To_Face(const T_GRID& grid,const int axis,const TV_INT& face_index,const T_ARRAYS_VECTOR& u_cell) const // this never needs to set starting points doesn't use velocities
    {FACE_ITERATOR iterator(grid,axis,face_index);
    TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face_index,cell1,cell2);
    if(body_list.Occupied_Face_Center(iterator)){
        T cell1_value=default_replacement_value,cell2_value=default_replacement_value;
        if(body_list.Cell_Center_Visible_From_Face(cell1,axis,face_index)) cell1_value=u_cell(cell1)[axis];
        if(body_list.Cell_Center_Visible_From_Face(cell2,axis,face_index)) cell2_value=u_cell(cell2)[axis];
        return (T).5*(cell1_value+cell2_value);}
    else return (T).5*(u_cell(cell1)[axis]+u_cell(cell2)[axis]);}

    TV Face_To_Face_Vector(const T_GRID& grid,const int side,const int axis,const TV_INT& face_index,const T_FACE_LOOKUP& u_face) const
    {FACE_ITERATOR iterator(grid,axis,face_index);
    const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(side,axis,face_index);
    return AVERAGING_UNIFORM<T_GRID,T_FACE_LOOKUP>::Average_Face_To_Face_Vector_Helper(grid,iterator,lookup);}

//#####################################################################
};
}
#endif
