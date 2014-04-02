//#####################################################################
// Copyright 2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AVERAGING_COLLIDABLE_DYADIC
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __AVERAGING_COLLIDABLE_DYADIC__
#define __AVERAGING_COLLIDABLE_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/AVERAGING_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP>
class AVERAGING_COLLIDABLE_DYADIC
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL T_CELL;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
public:
    typedef T_FACE_LOOKUP FACE_LOOKUP;

    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    T default_replacement_value;

    AVERAGING_COLLIDABLE_DYADIC(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const T default_replacement_value_input)
        :body_list(body_list_input),default_replacement_value(default_replacement_value_input)
    {}

    ~AVERAGING_COLLIDABLE_DYADIC()
    {}

    TV Face_To_Cell_Vector(const T_GRID& grid,const T_CELL* cell,const T_FACE_LOOKUP& u_face) const
    {TV value;const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell->Cell());lookup.Set_Reference_Point(cell->Center());
    for(int axis=1;axis<=T_GRID::dimension;axis++) value[axis]=(T).5*(lookup(cell->Face(2*axis-2))+lookup(cell->Face(2*axis-1)));
    return value;}

    T Cell_To_Face(const T_GRID& grid,const int face_index,const ARRAY<T>& u_cell,const ARRAY<TV>& u_node) const
    {FACE_ITERATOR iterator(grid,face_index);T_CELL *cell1=iterator.Deepest_Cell(),*cell2=iterator.Other_Cell();
    if(cell1->Depth_Of_This_Cell()==cell2->Depth_Of_This_Cell()){
        if(body_list.Occupied_Face_Center(iterator)){
            T cell1_value=default_replacement_value,cell2_value=default_replacement_value;
            if(body_list.Cell_Center_Visible_From_Face(cell1,iterator.Axis(),iterator.Face_Index_In_Cell())) cell1_value=u_cell(cell1->Cell());
            if(body_list.Cell_Center_Visible_From_Face(cell2,iterator.Axis(),iterator.Face_Index_In_Cell())) cell2_value=u_cell(cell2->Cell());
            return (T).5*(cell1_value+cell2_value);}
        else return (T).5*(u_cell(cell1->Cell())+u_cell(cell2->Cell()));}
    else{
        T sum=0;for(int node=0;node<T_GRID::number_of_nodes_per_face;node++) sum+=u_node(iterator.Face_Node(node));
        return T_GRID::one_over_number_of_nodes_per_face*sum;}}

    TV Face_To_Face_Vector(const T_GRID& grid,const int face_index,const T_FACE_LOOKUP& u_face) const
    {FACE_ITERATOR iterator(grid,face_index);assert(iterator.Other_Cell()->Depth_Of_This_Cell()==grid.maximum_depth);
    const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(face_index);lookup.Set_Reference_Point(iterator.Location());
    return AVERAGING_DYADIC<T_GRID,T_FACE_LOOKUP>::Average_Face_To_Face_Vector_Helper(grid,iterator,lookup);}

//#####################################################################
};
}
#endif
#endif
