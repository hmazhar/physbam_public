//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_UNIFORM
//##################################################################### 
#ifndef __VORTICITY_UNIFORM__
#define __VORTICITY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV>
class VORTICITY_UNIFORM
{
    typedef typename TV::SCALAR T;typedef typename REBIND<TV,int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename TV::SPIN T_SPIN;typedef typename T_ARRAYS_SCALAR::template REBIND<T_SPIN>::TYPE T_ARRAYS_SPIN;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;

private:
    template<class T_FACE_LOOKUP_LOOKUP>
    static T Cell_Centered_Face_Velocity_Difference(const GRID<TV>& grid,const int axis,const int difference_axis,const T_FACE_LOOKUP_LOOKUP& lookup,const TV_INT& cell_index)
    {T one_over_four_dh=(T).25*grid.one_over_dX[difference_axis];TV_INT difference_offset=TV_INT::Axis_Vector(difference_axis);
    TV_INT face_1=grid.First_Face_Index_In_Cell(axis,cell_index),face_2=grid.Second_Face_Index_In_Cell(axis,cell_index);
    return one_over_four_dh*(lookup(axis,face_1+difference_offset)+lookup(axis,face_2+difference_offset)-lookup(axis,face_1-difference_offset)-lookup(axis,face_2-difference_offset));}

public:
    template<class T_FACE_LOOKUP_LOOKUP>
    static T_SPIN Vorticity(const GRID<VECTOR<T,1> >& grid,const T_FACE_LOOKUP_LOOKUP& lookup,const TV_INT& index)
    {return T_SPIN();}

    template<class T_FACE_LOOKUP_LOOKUP>
    static T_SPIN Vorticity(const GRID<VECTOR<T,2> >& grid,const T_FACE_LOOKUP_LOOKUP& lookup,const TV_INT& index)
    {return VECTOR<T,1>(Cell_Centered_Face_Velocity_Difference(grid,2,1,lookup,index)-Cell_Centered_Face_Velocity_Difference(grid,1,2,lookup,index));}
    template<class T_FACE_LOOKUP_LOOKUP>
    static T_SPIN Vorticity(const GRID<VECTOR<T,3> >& grid,const T_FACE_LOOKUP_LOOKUP& lookup,const TV_INT& index)
    {return VECTOR<T,3>(
        Cell_Centered_Face_Velocity_Difference(grid,3,2,lookup,index)-Cell_Centered_Face_Velocity_Difference(grid,2,3,lookup,index),
        Cell_Centered_Face_Velocity_Difference(grid,1,3,lookup,index)-Cell_Centered_Face_Velocity_Difference(grid,3,1,lookup,index),
        Cell_Centered_Face_Velocity_Difference(grid,2,1,lookup,index)-Cell_Centered_Face_Velocity_Difference(grid,1,2,lookup,index));}

    template<class T_FACE_LOOKUP_2>
    static void Vorticity(const GRID<TV>& grid,const T_FACE_LOOKUP_2& face_velocities_lookup,T_ARRAYS_SPIN& vorticity,T_ARRAYS_SCALAR& vorticity_magnitude)
    {for(CELL_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        const typename T_FACE_LOOKUP_2::LOOKUP& lookup=face_velocities_lookup.Starting_Point_Cell(iterator.Cell_Index());lookup.Set_Reference_Point(iterator.Location());
        vorticity(index)=Vorticity(grid,lookup,index);vorticity_magnitude(index)=vorticity(index).Magnitude();}}
};
}
#endif
