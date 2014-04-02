//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAYS
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ARRAYS__
#define __READ_WRITE_ARRAYS__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T,int dimension>
class Read_Write<ARRAYS_ND_BASE<VECTOR<T,dimension> >,RW>
{
    typedef VECTOR<T,dimension> TV;typedef VECTOR<int,dimension> TV_INT;
    typedef ARRAYS_ND_BASE<TV> T_ARRAYS;
public:
    static void Read(std::istream& input,T_ARRAYS& object)
    {Read_With_Length(input,1,object);}

    static void Write(std::ostream& output,const T_ARRAYS& object)
    {Write_With_Length(output,1,object);}

protected:
    static void Read_With_Length(std::istream& input,const int length2,T_ARRAYS& object)
    {int read_length;Read_Binary<RW>(input,read_length,object.domain);
    if(read_length!=length2) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Read length %d not equal to %d",read_length,length2));
    if(object.counts.Min()<0) throw READ_ERROR("Invalid negative array size");
    object.counts=object.domain.Edge_Lengths()+1;
    int size=object.counts.Product();
    if(size!=object.array.Size()){delete[] object.array.Get_Array_Pointer();ARRAY_VIEW<T> new_array(size,new T[size]);object.array.Exchange(new_array);}
    Read_Binary_Array<RW>(input,object.array.Get_Array_Pointer(),object.array.Size());object.Calculate_Acceleration_Constants();}

    static void Write_With_Length(std::ostream& output,const int length2,const T_ARRAYS& object)
    {Write_Binary<RW>(output,length2,object.domain);Write_Binary_Array<RW>(output,object.array.Get_Array_Pointer(),object.array.Size());}

};

template<class RW,class T,int dimension>
class Read_Write<ARRAY<T,VECTOR<int,dimension> >,RW>:public Read_Write<ARRAYS_ND_BASE<VECTOR<T,dimension> >,RW>
{};

template<class RW,class T,int dimension>
class Read_Write<ARRAY_VIEW<T,VECTOR<int,dimension> >,RW>:public Read_Write<ARRAYS_ND_BASE<VECTOR<T,dimension> >,RW>
{
public:
    static void Read(std::istream& input,ARRAY_VIEW<T,VECTOR<int,dimension> >& object)
    {PHYSBAM_NOT_IMPLEMENTED();}
};

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAYS_ND_BASE<VECTOR<T,1> >& a)
{for(int i=a.domain.min_corner.x;i<=a.domain.max_corner.x;i++) output_stream<<a(i)<<" ";output_stream<<std::endl;return output_stream;}

template<class T> inline std::ostream& operator<<(std::ostream& output,const ARRAYS_ND_BASE<VECTOR<T,2> >& a)
{for(int i=a.domain.min_corner.x;i<=a.domain.max_corner.x;i++){for(int j=a.domain.min_corner.y;j<=a.domain.max_corner.y;j++) output<<a(i,j)<<" ";output<<std::endl;}return output;}

template<class T> inline std::ostream& operator<<(std::ostream& output,const ARRAYS_ND_BASE<VECTOR<T,3> >& a)
{for(int i=a.domain.min_corner.x;i<=a.domain.max_corner.x;i++){for(int j=a.domain.min_corner.y;j<=a.domain.max_corner.y;j++){for(int ij=a.domain.min_corner.z;ij<=a.domain.max_corner.z;ij++)output<<a(i,j,ij)<<" ";output<<std::endl;}
    output<<"------------------------------------------"<<std::endl;}
return output;}
}
#endif
#endif
