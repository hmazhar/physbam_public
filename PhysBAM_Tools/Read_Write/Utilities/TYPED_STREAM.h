//#####################################################################
// Copyright 2006, Geoffrey Irving, Eftychios Sifakis, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TYPED_ISTREAM and TYPED_OSTREAM
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __TYPED_STREAM__
#define __TYPED_STREAM__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <iostream>
namespace PhysBAM{

//#####################################################################
// Class STREAM_TYPE
//#####################################################################
class STREAM_TYPE
{
public:
    const bool use_doubles; // otherwise use floats

    explicit STREAM_TYPE(const float)
        :use_doubles(false)
    {}

    explicit STREAM_TYPE(const double)
        :use_doubles(true)
    {}
};
//#####################################################################
// Class TYPED_ISTREAM
//#####################################################################
class TYPED_ISTREAM
{
public:
    std::istream& stream;
    const STREAM_TYPE type;

    TYPED_ISTREAM(std::istream& stream_input,const STREAM_TYPE type_input)
        :stream(stream_input),type(type_input)
    {}

    template<class T> void Read(T& d)
    {if(!type.use_doubles) Read_Write<T,float>::Read(stream,d);
    else Read_Write<T,double>::Read(stream,d);}
};
//#####################################################################
// Class TYPED_OSTREAM
//#####################################################################
class TYPED_OSTREAM
{
public:
    std::ostream& stream;
    const STREAM_TYPE type;

    TYPED_OSTREAM(std::ostream& stream_input,const STREAM_TYPE type_input)
        :stream(stream_input),type(type_input)
    {}

    template<class T> void Write(const T& d)
    {if(!type.use_doubles) Read_Write<T,float>::Write(stream,d);
    else Read_Write<T,double>::Write(stream,d);}
};
//#####################################################################
// Detect whether a type has Read/Write taking typed streams
//#####################################################################
template<class T,class HAS> struct HAS_TYPED_READ{enum {value=false};};
template<class T> struct HAS_TYPED_READ<T,typename FIRST<void,typename T::HAS_TYPED_READ_WRITE>::TYPE>{enum {value=true};};
template<class T,class HAS> struct HAS_TYPED_WRITE{enum {value=false};};
template<class T> struct HAS_TYPED_WRITE<T,typename FIRST<void,typename T::HAS_TYPED_READ_WRITE>::TYPE>{enum {value=true};};
//#####################################################################
}
#endif
#endif
