//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FACE_ARRAYS
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FACE_ARRAYS__
#define __READ_WRITE_FACE_ARRAYS__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<ARRAY<T,FACE_INDEX<d> >,RW>
{
    typedef VECTOR<int,d> TV_INT;
    typedef ARRAY_VIEW<T,TV_INT> T_ARRAY_VIEW;
public:
    static void Read(std::istream& input,ARRAY<T,FACE_INDEX<d> >& object)
    {object.Clean_Memory();Read_Binary<RW>(input,object.domain_indices);Read_Binary<RW>(input,object.buffer_size);
    if(object.buffer_size<0) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid negative array size %d",object.buffer_size));
    if(!object.buffer_size) return;
    object.base_pointer=new T[object.buffer_size];
    Read_Binary_Array<RW>(input,object.base_pointer,object.buffer_size);
    T* p_start=object.base_pointer;
    for(int i=1;i<=d;i++){
        RANGE<TV_INT> domain;
        Read_Binary<RW>(input,domain);
        T_ARRAY_VIEW array_new(domain,p_start);
        T_ARRAY_VIEW::Exchange_Arrays(array_new,object.data(i));
        p_start+=(domain.Edge_Lengths()+1).Product();}}

    static void Write(std::ostream& output,const ARRAY<T,FACE_INDEX<d> >& object)
    {Write_Binary<RW>(output,object.domain_indices);Write_Binary<RW>(output,object.buffer_size);Write_Binary_Array<RW>(output,object.base_pointer,object.buffer_size);for(int i=1;i<=d;i++) Write_Binary<RW>(output,object.data(i).domain);}
};

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAY<T,FACE_INDEX<1> >& a)
{for(int i=a.domain_indices.min_corner.x;i<=a.domain_indices.max_corner.x+1;i++) output_stream<<a.Component(1)(i)<<" ";output_stream<<std::endl;return output_stream;}

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAY<T,FACE_INDEX<2> >& a)
{PHYSBAM_NOT_IMPLEMENTED();}

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAY<T,FACE_INDEX<3> >& a)
{PHYSBAM_NOT_IMPLEMENTED();}

}
#endif
#endif
