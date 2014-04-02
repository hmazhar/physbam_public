//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// class FACE_LOOKUP_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __FACE_LOOKUP_RLE__
#define __FACE_LOOKUP_RLE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

template<class T_GRID>
class FACE_LOOKUP_RLE
{
    typedef typename T_GRID::SCALAR T;
public:
    typedef T ELEMENT;

    const ARRAY<T>& V_face;

    FACE_LOOKUP_RLE(const ARRAY<T>& V_face_input)
        :V_face(V_face_input)
    {}

    typedef FACE_LOOKUP_RLE LOOKUP;

    const ARRAY<T>& Raw_Data() const
    {return V_face;}

    template<class T_FACE>
    const LOOKUP& Starting_Point_Face(const T_FACE& face) const
    {return *this;}

    const LOOKUP& Starting_Point_Cell(const int cell) const
    {return *this;}

    T operator()(const int axis,const int face) const
    {return V_face(face);}

//#####################################################################
};
}
#endif
#endif
