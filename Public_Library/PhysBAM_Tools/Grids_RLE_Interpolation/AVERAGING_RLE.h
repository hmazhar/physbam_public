//#####################################################################
// Copyright 2005, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __AVERAGING_RLE__
#define __AVERAGING_RLE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_RLE<T_GRID>
class AVERAGING_RLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::BLOCK T_BLOCK;
public:
    typedef T_FACE_LOOKUP FACE_LOOKUP;

    AVERAGING_RLE()
    {}

    ~AVERAGING_RLE()
    {}

    TV Face_To_Cell_Vector(const CELL_ITERATOR& cell,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell);
    TV value;for(int axis=1;axis<=T_GRID::dimension;axis++)
        value[axis]=(T).5*(lookup(axis,cell.First_Face_Index(axis))+lookup(axis,cell.Second_Face_Index(axis)));
    return value;}

    TV Face_To_Cell_Vector(const T_BLOCK& block,const CELL_ITERATOR& cell,const T_FACE_LOOKUP& u_face) const
    {return Face_To_Cell_Vector(cell,u_face);}

    // Note this never has to be replaced because it is phi averaged to the faces, not the velocities
    template<class T_FACE>
    T Cell_To_Face(const int axis,const T_FACE& face,const ARRAY<T>& u_cell) const
    {return (T).5*(u_cell(face.First_Cell_Index())+u_cell(face.Second_Cell_Index()));}

    template<class T_FACE>
    TV Face_To_Face_Vector(const T_FACE& face,const T_FACE_LOOKUP& u_face) const
    {return Average_Face_To_Face_Vector(face,u_face);}

    template<class T_FACE>
    TV Face_To_Face_Vector(const T_BLOCK& block,const T_FACE& face,const T_FACE_LOOKUP& u_face) const
    {return Face_To_Face_Vector(face,u_face);}

    template<class T_FACE_ITERATOR,class T_FACE_LOOKUP_LOOKUP>
    static typename T_GRID::VECTOR_T Average_Face_To_Face_Vector(const T_FACE_ITERATOR& iterator,T_FACE_LOOKUP_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(iterator);
    return Average_Face_To_Face_Vector_Helper(iterator,lookup,typename T_GRID::GRID_TAG());}

    template<class T_FACE_ITERATOR,class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,1> Average_Face_To_Face_Vector_Helper(const T_FACE_ITERATOR& iterator,const T_FACE_LOOKUP_LOOKUP& u_face,RLE_TAG<VECTOR<T,1> >)
    {PHYSBAM_NOT_IMPLEMENTED();}

    template<class T_FACE_ITERATOR,class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,2> Average_Face_To_Face_Vector_Helper(const T_FACE_ITERATOR& iterator,const T_FACE_LOOKUP_LOOKUP& u_face,RLE_TAG<VECTOR<T,2> >)
    {int axis=iterator.Axis();typename RLE_GRID_2D<T>::INDEX face=iterator.Face();
    TV value;value[axis]=u_face(axis,face);
    int other_axis=3-axis;
    value[other_axis]=(T).25*(u_face(other_axis,iterator.cell1.First_Face_Index(other_axis))+u_face(other_axis,iterator.cell1.Second_Face_Index(other_axis))+
        u_face(other_axis,iterator.cell2.First_Face_Index(other_axis))+u_face(other_axis,iterator.cell2.Second_Face_Index(other_axis)));
    return value;}

    template<class T_FACE_ITERATOR,class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,3> Average_Face_To_Face_Vector_Helper(const T_FACE_ITERATOR& iterator,const T_FACE_LOOKUP_LOOKUP& u_face,RLE_TAG<VECTOR<T,3> >)
    {static const int axis_to_other_axis[3][2]={{2,3},{1,3},{1,2}};
    int axis=iterator.Axis();typename RLE_GRID_3D<T>::INDEX face=iterator.Face();
    VECTOR<T,3> value;value[axis]=u_face(axis,face);
    for(int i=1;i<=T_GRID::dimension-1;i++){
        int other_axis=axis_to_other_axis[axis-1][i-1];
        value[other_axis]=(T).25*(u_face(other_axis,iterator.cell1.First_Face_Index(other_axis))+u_face(other_axis,iterator.cell1.Second_Face_Index(other_axis))+
                            u_face(other_axis,iterator.cell2.First_Face_Index(other_axis))+u_face(other_axis,iterator.cell2.Second_Face_Index(other_axis)));}
    return value;}

//#####################################################################
};
}
#endif
#endif
