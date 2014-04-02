//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_FACE_INDEX
//#####################################################################
#ifndef __DYADIC_FACE_INDEX__
#define __DYADIC_FACE_INDEX__

namespace PhysBAM{

template<int d>
class DYADIC_FACE_INDEX
{
    typedef int INDEX;
public:
    DYADIC_FACE_INDEX()
        :axis(0)
    {}

    DYADIC_FACE_INDEX(int axis_input,const INDEX& index_input,const INDEX& first_cell_index,const INDEX& second_cell_index)
        :axis(axis_input),index(index_input),first_cell_index(first_cell_index),second_cell_index(second_cell_index)
    {}

    int axis;
    INDEX index,first_cell_index,second_cell_index;

    bool operator==(const DYADIC_FACE_INDEX& fi) const
    {return axis==fi.axis && index==fi.index;}

    INDEX First_Cell_Index() const
    {return first_cell_index;}

    INDEX Second_Cell_Index() const
    {return second_cell_index;}

    // INDEX Cell_Index(int j) const
    // {INDEX i(index);i(axis)+=j-2;return i;}
};
}
// #include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_DYADIC_FACE_INDEX.h>
#endif
