//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_UNIFORM
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class RW,class T_GRID> void Read_Write<LEVELSET_UNIFORM<T_GRID>,RW>::
Read(std::istream& input,LEVELSET_UNIFORM<T_GRID>& object)
{
    Read_Binary<RW>(input,object.grid,object.phi);
}
//#####################################################################
// Function Read
//#####################################################################
template<class RW,class T_GRID> void Read_Write<LEVELSET_UNIFORM<T_GRID>,RW>::
Write(std::ostream& output,const LEVELSET_UNIFORM<T_GRID>& object)
{
    Write_Binary<RW>(output,object.grid,object.phi);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,RW) \
    template class Read_Write<LEVELSET_UNIFORM<GRID<VECTOR<T,1> > >,RW>; \
    template class Read_Write<LEVELSET_UNIFORM<GRID<VECTOR<T,2> > >,RW>; \
    template class Read_Write<LEVELSET_UNIFORM<GRID<VECTOR<T,3> > >,RW>;
INSTANTIATION_HELPER(float,float);
INSTANTIATION_HELPER(float,double);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,float);
INSTANTIATION_HELPER(double,double);
#endif

#endif
