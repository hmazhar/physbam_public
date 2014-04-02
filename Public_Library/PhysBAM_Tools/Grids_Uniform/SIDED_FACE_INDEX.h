//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIDED_FACE_INDEX
//#####################################################################
#ifndef __SIDED_FACE_INDEX__
#define __SIDED_FACE_INDEX__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<int d>
class SIDED_FACE_INDEX
{
    typedef VECTOR<int,d> TV_INT;
public:
    SIDED_FACE_INDEX()
        :side(0),axis(0)
    {}

    SIDED_FACE_INDEX(int side_input,int axis_input,const TV_INT& index_input)
        :side(side_input),axis(axis_input),index(index_input)
    {}

    SIDED_FACE_INDEX(int side_input,const FACE_INDEX<d>& face_index)
        :side(side_input),axis(face_index.axis),index(face_index.index)
    {}

    int side;
    int axis;
    TV_INT index;

    bool operator==(const SIDED_FACE_INDEX& fi) const
    {return side==fi.side && axis==fi.axis && index==fi.index;}

    FACE_INDEX<d> Face_Index() const
    {return FACE_INDEX<d>(axis,index);}

    TV_INT First_Cell_Index() const
    {TV_INT i(index);i(axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return index;}

    TV_INT Real_Cell_Index() const
    {return side==1?index-TV_INT::Axis_Vector(axis):index;}

    TV_INT Ghost_Cell_Index() const
    {return side==2?index-TV_INT::Axis_Vector(axis):index;}
};
}
#endif
