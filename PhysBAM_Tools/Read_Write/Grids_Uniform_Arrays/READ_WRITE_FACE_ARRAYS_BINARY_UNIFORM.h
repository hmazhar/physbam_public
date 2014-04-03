//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FACE_ARRAYS_BINARY_UNIFORM
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FACE_ARRAYS_BINARY_UNIFORM__
#define __READ_WRITE_FACE_ARRAYS_BINARY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<ARRAY<T,SIDED_FACE_INDEX<d> >,RW>
{
    typedef VECTOR<int,d> TV_INT;
    typedef ARRAY_VIEW<T,TV_INT> T_ARRAY_VIEW;
public:
    static void Read(std::istream& input,ARRAY<T,SIDED_FACE_INDEX<d> >& object)
    {Read_Binary<RW>(input,dynamic_cast<ARRAY<T,FACE_INDEX<d> >&>(object),object.u2);object.Initialize();}

    static void Write(std::ostream& output,const ARRAY<T,SIDED_FACE_INDEX<d> >& object)
    {Write_Binary<RW>(output,dynamic_cast<const ARRAY<T,FACE_INDEX<d> >&>(object),object.u2);}
};
}
#endif
#endif
