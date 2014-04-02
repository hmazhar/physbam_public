//#####################################################################
// Copyright 2005, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// class FACE_LOOKUP_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __FACE_LOOKUP_DYADIC__
#define __FACE_LOOKUP_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

template<class T_GRID>
class FACE_LOOKUP_DYADIC
{
    typedef typename T_GRID::SCALAR T;
public:
    typedef T ELEMENT;

    const ARRAY<T>& V_face;

    FACE_LOOKUP_DYADIC(const ARRAY<T>& V_face_input)
        :V_face(V_face_input)
    {}

    typedef FACE_LOOKUP_DYADIC LOOKUP;

    const ARRAY<T>& Raw_Data() const
    {return V_face;}

    const LOOKUP& Starting_Point_Face(const int face) const
    {return *this;}

    const LOOKUP& Starting_Point_Cell(const int cell) const
    {return *this;}

    T operator()(const int face) const
    {return V_face(face);}

//#####################################################################
};
}
#endif
#endif
