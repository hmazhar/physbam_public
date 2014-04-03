//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_LOOKUP_UNIFORM
//#####################################################################
#ifndef __FACE_LOOKUP_UNIFORM__
#define __FACE_LOOKUP_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class FACE_LOOKUP_UNIFORM
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::INDEX T_INDEX;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;
public:
    typedef T ELEMENT;
        
    const T_FACE_ARRAYS& V_face;

    FACE_LOOKUP_UNIFORM(const T_FACE_ARRAYS& V_face_input)
        :V_face(V_face_input)
    {}

    const T_FACE_ARRAYS& Raw_Data() const
    {return V_face;}

    int Number_Of_Ghost_Cells() const
    {return V_face.Number_Of_Ghost_Cells();}
    
    typedef FACE_LOOKUP_UNIFORM LOOKUP;

    const LOOKUP& Starting_Point_Face(const int axis,const T_INDEX& face) const
    {return *this;}
    
    const LOOKUP& Starting_Point_Cell(const T_INDEX& cell) const
    {return *this;}

    void Set_Reference_Point(const TV& reference_point) const
    {}

    T operator()(const int axis,const T_INDEX& face) const
    {return V_face.Component(axis)(face);}

    T operator()(const FACE_INDEX<TV::dimension>& face) const
    {return V_face(face);}

//#####################################################################
};
}
#endif
