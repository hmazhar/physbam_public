//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ROTATION
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ROTATION__
#define __READ_WRITE_ROTATION__

#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_COMPLEX.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_QUATERNION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<ROTATION<VECTOR<T,1> >,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<ROTATION<VECTOR<T,1> >,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,ROTATION<VECTOR<T,1> >& object)
    {}

    static void Write(std::ostream& output,const ROTATION<VECTOR<T,1> >& object)
    {}
};

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const ROTATION<VECTOR<T,1> >& r)
{return output_stream;}

template<class T>
inline std::istream& operator>>(std::istream& input_stream,ROTATION<VECTOR<T,1> >& r)
{return input_stream;}

template<class RW,class T>
class Read_Write<ROTATION<VECTOR<T,2> >,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<ROTATION<VECTOR<T,2> >,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,ROTATION<VECTOR<T,2> >& object)
    {Read_Binary<RW>(input,object.c);}

    static void Write(std::ostream& output,const ROTATION<VECTOR<T,2> >& object)
    {Write_Binary<RW>(output,object.c);}
};

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const ROTATION<VECTOR<T,2> >& r)
{return output_stream<<r.Complex();}

template<class T>
inline std::istream& operator>>(std::istream& input_stream,ROTATION<VECTOR<T,2> >& r)
{COMPLEX<T> c;input_stream>>c;r=ROTATION<VECTOR<T,3> >::From_Complex(c);return input_stream;}


template<class RW,class T>
class Read_Write<ROTATION<VECTOR<T,3> >,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<ROTATION<VECTOR<T,3> >,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,ROTATION<VECTOR<T,3> >& object)
    {Read_Binary<RW>(input,object.q);if(!IS_SAME<T,int>::value && !object.Is_Normalized()) PHYSBAM_FATAL_ERROR("Read nonnormalized rotation");}

    static void Write(std::ostream& output,const ROTATION<VECTOR<T,3> >& object)
    {Write_Binary<RW>(output,object.q);}
};

template<class T>
inline std::ostream& operator<<(std::ostream& output,const ROTATION<VECTOR<T,3> >& r)
{return output<<r.Quaternion();}

template<class T>
inline std::istream& operator>>(std::istream& input,ROTATION<VECTOR<T,3> >& r)
{QUATERNION<T> q;input>>q;r=ROTATION<VECTOR<T,3> >::From_Quaternion(q);return input;}
}
#endif
#endif
