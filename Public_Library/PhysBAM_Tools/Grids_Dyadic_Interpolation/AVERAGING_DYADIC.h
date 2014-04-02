//#####################################################################
// Copyright 2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AVERAGING_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __AVERAGING_DYADIC__
#define __AVERAGING_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID>
class AVERAGING_DYADIC
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::CELL T_CELL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
public:
    typedef T_FACE_LOOKUP FACE_LOOKUP;

    AVERAGING_DYADIC()
    {}

    ~AVERAGING_DYADIC()
    {}

    TV Face_To_Cell_Vector(const T_GRID& grid,const T_CELL* cell,const T_FACE_LOOKUP& u_face) const
    {TV value;const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell->Cell());
    for(int axis=1;axis<=T_GRID::dimension;axis++) value[axis]=(T).5*(lookup(cell->Face(2*axis-2))+lookup(cell->Face(2*axis-1)));
    return value;}

    // Note this never has to be replaced because it is phi averaged to the faces, not the velocities
    T Cell_To_Face(const T_GRID& grid,const int face_index,const ARRAY<T>& u_cell,const ARRAY<T>& u_node) const
    {FACE_ITERATOR iterator(grid,face_index);T_CELL *cell1=iterator.Deepest_Cell(),*cell2=iterator.Other_Cell();
    if(cell1->Depth_Of_This_Cell()==cell2->Depth_Of_This_Cell()) return (T).5*(u_cell(cell1->Cell())+u_cell(cell2->Cell()));
    else{
        T sum=0;for(int node=0;node<T_GRID::number_of_nodes_per_face;node++)sum+=u_node(iterator.Face_Node(node));
        return T_GRID::one_over_number_of_nodes_per_face*sum;}}

    TV Face_To_Face_Vector(const T_GRID& grid,const int face_index,const T_FACE_LOOKUP& u_face,const ARRAY<TV>& u_node) const
    {FACE_ITERATOR iterator(grid,face_index);
    if(iterator.Other_Cell()->Depth_Of_This_Cell()==grid.maximum_depth) return Average_Face_To_Face_Vector(grid,iterator,u_face);
    return Face_To_Face_Vector_Helper(grid,iterator.Deepest_Cell(),u_face.V_face,u_node,iterator.Axis(),iterator.Face_Index_In_Cell());}

private:
    static TV Face_To_Face_Vector_Helper(const BINTREE_GRID<T>& grid,const BINTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const int axis,const int face_index_in_cell)
    {assert(0<=face_index_in_cell&&face_index_in_cell<2);
     return TV(u_face(cell->Face(face_index_in_cell)));
    }

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT    
    static TV Face_To_Face_Vector_Helper(const QUADTREE_GRID<T>& grid,const QUADTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const int axis,const int face_index_in_cell)
    {assert(0<=face_index_in_cell&&face_index_in_cell<4);static const bool cell_position[6]={1,0,1,0,1,0};
    int node1=cell->Node(MAP_QUADTREE_MESH<T>::nodes_by_axis[axis][cell_position[face_index_in_cell]][0]),
        node2=cell->Node(MAP_QUADTREE_MESH<T>::nodes_by_axis[axis][cell_position[face_index_in_cell]][1]);
    switch(axis){
        case 0:return TV(u_face(cell->Face(face_index_in_cell)),(T).5*(u_node(node1).y+u_node(node2).y));
        default:assert(axis==1);return TV((T).5*(u_node(node1).x+u_node(node2).x),u_face(cell->Face(face_index_in_cell)));}}

    static TV Face_To_Face_Vector_Helper(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const int axis,const int face_index_in_cell)
    {assert(0<=face_index_in_cell&&face_index_in_cell<6);static const bool cell_position[6]={1,0,1,0,1,0};
    int node1=cell->Node(MAP_OCTREE_MESH<T>::nodes_by_axis[axis][cell_position[face_index_in_cell]][0]),node2=cell->Node(MAP_OCTREE_MESH<T>::nodes_by_axis[axis][cell_position[face_index_in_cell]][1]),
        node3=cell->Node(MAP_OCTREE_MESH<T>::nodes_by_axis[axis][cell_position[face_index_in_cell]][2]),node4=cell->Node(MAP_OCTREE_MESH<T>::nodes_by_axis[axis][cell_position[face_index_in_cell]][3]);
    switch(axis){
        case 0:return TV(u_face(cell->Face(face_index_in_cell)),(T).25*(u_node(node1).y+u_node(node2).y+u_node(node3).y+u_node(node4).y),
            (T).25*(u_node(node1).z+u_node(node2).z+u_node(node3).z+u_node(node4).z));
        case 1:return TV((T).25*(u_node(node1).x+u_node(node2).x+u_node(node3).x+u_node(node4).x),u_face(cell->Face(face_index_in_cell)),
            (T).25*(u_node(node1).z+u_node(node2).z+u_node(node3).z+u_node(node4).z));
        default:assert(axis==2);return TV((T).25*(u_node(node1).x+u_node(node2).x+u_node(node3).x+u_node(node4).x),
            (T).25*(u_node(node1).y+u_node(node2).y+u_node(node3).y+u_node(node4).y),u_face(cell->Face(face_index_in_cell)));}}
#endif

public:
    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,1> Average_Face_To_Face_Vector(const BINTREE_GRID<T>& grid,const DYADIC_GRID_ITERATOR_FACE<BINTREE_GRID<T> >& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {return VECTOR<T,1>(u_face(iterator.Face_Index()));}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    template<class T_FACE_ITERATOR>
    static VECTOR<T,2> Average_Face_To_Face_Vector(const QUADTREE_GRID<T>& grid,const T_FACE_ITERATOR& iterator,const T_FACE_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(iterator.Face_Index());
    return Average_Face_To_Face_Vector_Helper(grid,iterator,lookup);}

    template<class T_FACE_ITERATOR>
    static VECTOR<T,3> Average_Face_To_Face_Vector(const OCTREE_GRID<T>& grid,const T_FACE_ITERATOR& iterator,const T_FACE_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(iterator.Face_Index());
    return Average_Face_To_Face_Vector_Helper(grid,iterator,lookup);}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,2> Average_Face_To_Face_Vector_Helper(const QUADTREE_GRID<T>& grid,const DYADIC_GRID_ITERATOR_FACE<QUADTREE_GRID<T> >& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {int axis=iterator.Axis();typename QUADTREE_GRID<T>::INDEX face=iterator.Face_Index(),cell1,cell2;iterator.Unordered_Cell_Indices_Touching_Face(cell1,cell2);
    TV value;value[axis+1]=u_face(face);
    int other_axis=1-axis;
    value[other_axis+1]=(T).25*(u_face(grid.First_Face_Index_In_Cell(other_axis,cell1))+u_face(grid.Second_Face_Index_In_Cell(other_axis,cell1))+
                        u_face(grid.First_Face_Index_In_Cell(other_axis,cell2))+u_face(grid.Second_Face_Index_In_Cell(other_axis,cell2)));
    return value;}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,3> Average_Face_To_Face_Vector_Helper(const OCTREE_GRID<T>& grid,const DYADIC_GRID_ITERATOR_FACE<OCTREE_GRID<T> >& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {static const int axis_to_other_axis[3][2]={{1,2},{0,2},{0,1}};
    int axis=iterator.Axis();typename OCTREE_GRID<T>::INDEX face=iterator.Face_Index(),cell1,cell2;iterator.Unordered_Cell_Indices_Touching_Face(cell1,cell2);
    VECTOR<T,3> value;value[axis+1]=u_face(face);
    for(int i=1;i<=T_GRID::dimension-1;i++){
        int other_axis=axis_to_other_axis[axis][i-1];
        value[other_axis+1]=(T).25*(u_face(grid.First_Face_Index_In_Cell(other_axis,cell1))+u_face(grid.Second_Face_Index_In_Cell(other_axis,cell1))+
                            u_face(grid.First_Face_Index_In_Cell(other_axis,cell2))+u_face(grid.Second_Face_Index_In_Cell(other_axis,cell2)));}
    return value;}
#endif
//#####################################################################
};
}
#endif
#endif
