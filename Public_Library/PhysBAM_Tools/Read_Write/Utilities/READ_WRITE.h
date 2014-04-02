//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class Read_Write
//#####################################################################
// Functions for reading and writing which do the correct thing for objects, pointers, primitive types, etc. In general, use Read/Write_Binary (and Read/Write_Binary_Array) using T for the type
// of the object you're reading/writing and RW the underlying floating point scalar type (float/double).
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __Read_Write__
#define __Read_Write__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <cassert>
#include <iostream>
namespace PhysBAM{

#ifdef PHYSBAM_BIG_ENDIAN
static const bool big_endian=true;
#else
static const bool big_endian=false;
#endif

template<> struct PLATFORM_INDEPENDENT_SIZE<bool>{static const int value=1;};
template<> struct PLATFORM_INDEPENDENT_SIZE<char>{static const int value=1;};
template<> struct PLATFORM_INDEPENDENT_SIZE<unsigned char>{static const int value=1;};
template<> struct PLATFORM_INDEPENDENT_SIZE<short>{static const int value=2;};
template<> struct PLATFORM_INDEPENDENT_SIZE<unsigned short>{static const int value=2;};
template<> struct PLATFORM_INDEPENDENT_SIZE<int>{static const int value=4;};
template<> struct PLATFORM_INDEPENDENT_SIZE<unsigned int>{static const int value=4;};
template<> struct PLATFORM_INDEPENDENT_SIZE<float>{static const int value=4;};
template<> struct PLATFORM_INDEPENDENT_SIZE<double>{static const int value=8;};
template<class T> struct PLATFORM_INDEPENDENT_SIZE<T,typename ENABLE_IF<IS_ENUM<T>::value>::TYPE>{static const int value=4;};

template<class T> struct IS_PRIMITIVE_BINARY_IO_SAFE {static const bool value=!big_endian && (sizeof(T)==PLATFORM_INDEPENDENT_SIZE<T>::value);};

// classify types which can be written directly to files without conversion
template<class T,class RW,class ENABLER> struct IS_BINARY_IO_SAFE {static const bool value=false;};
template<class T,class RW> struct IS_BINARY_IO_SAFE<T,RW,typename ENABLE_IF<OR<IS_INTEGRAL<T>::value,IS_ENUM<T>::value>::value>::TYPE>:public IS_PRIMITIVE_BINARY_IO_SAFE<T>{};
template<> struct IS_BINARY_IO_SAFE<float,float>:public IS_PRIMITIVE_BINARY_IO_SAFE<float>{};
template<> struct IS_BINARY_IO_SAFE<double,double>:public IS_PRIMITIVE_BINARY_IO_SAFE<double>{};

//#####################################################################
// Function Swap_Endianity
//#####################################################################
template<class T>
inline void Swap_Endianity(T& x)
{assert(sizeof(T)<=8);
if(sizeof(T)>1) {T old=x;for(unsigned int k=1;k<=sizeof(T);k++) ((char*)&x)[k-1]=((char*)&old)[sizeof(T)-k];}}

//#####################################################################
// Read and Write for primitives
//#####################################################################
template<class T> inline void
Read_Primitive(std::istream& input,T& d)
{
    STATIC_ASSERT((sizeof(T)==PLATFORM_INDEPENDENT_SIZE<T>::value));
    input.read((char*)&d,sizeof(T));
    if(big_endian) Swap_Endianity(d); // convert to little endian if necessary
}

template<class T> inline void
Write_Primitive(std::ostream& output,const T& d)
{
    STATIC_ASSERT((sizeof(T)==PLATFORM_INDEPENDENT_SIZE<T>::value));
    if(big_endian){T d2=d;Swap_Endianity(d2);output.write((const char*)&d2,sizeof(T));} // convert to big endian if necessary
    else output.write((const char*)&d,sizeof(T));
}
//#####################################################################
// PhysBAM data formats assume sizeof(bool)==1 but the Mac apparently has bool==int with sizeof(bool)==4, so need to specialize these
//#####################################################################
#ifdef __APPLE__

template<class RW> struct IS_BINARY_IO_SAFE<bool,RW> {static const bool value=false;};

template<>
inline void Read_Primitive<bool>(std::istream& input,bool& d)
{char c;input.read(&c,1);d=(bool)c;}

template<>
inline void Write_Primitive<bool>(std::ostream& output,const bool& d)
{char c=(char)d;output.write(&c,1);}

#endif

//#####################################################################
// Read_Write for objects
//#####################################################################
template<class T,class RW,class ENABLER>
class Read_Write
{
public:
    // Making this "virtual" base class.  Nothing should use it directly.
    static void Read(std::istream& input,T& d)
    {STATIC_ASSERT((T)false);}

    static void Write(std::ostream& output,const T& d)
    {STATIC_ASSERT((T)false);}
};
//#####################################################################
// Read_Write for binary I/O safe types
//#####################################################################
template<class T,class RW> class Read_Write<T,RW,typename ENABLE_IF<IS_BINARY_IO_SAFE<T,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,T& d)
    {input.read(reinterpret_cast<char*>(&d),sizeof(T));}

    static void Write(std::ostream& output,const T& d)
    {output.write(reinterpret_cast<const char*>(&d),sizeof(T));}
};
//#####################################################################
// Read_Write for integral types (bool, char, int, etc.) and enums
//#####################################################################
template<class T,class RW> class Read_Write<T,RW,typename ENABLE_IF<(IS_INTEGRAL<T>::value || IS_ENUM<T>::value) && !IS_BINARY_IO_SAFE<T,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,T& d)
    {Read_Primitive(input,d);}

    static void Write(std::ostream& output,const T& d)
    {Write_Primitive(output,d);}
};
//#####################################################################
// Read_Write for float and double
//#####################################################################
template<class T,class RW> class Read_Write<T,RW,typename ENABLE_IF<(IS_SAME<T,float>::value || IS_SAME<T,double>::value) && !IS_BINARY_IO_SAFE<T,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,T& d)
    {RW tmp;Read_Primitive(input,tmp);d=(T)tmp;}

    static void Write(std::ostream& output,const T& d)
    {Write_Primitive(output,(RW)d);}
};
//#####################################################################
// Read_Write for pointers to data
//#####################################################################
template<class T,class RW>
class Read_Write<T*,RW>
{
public:
    static void Read(std::istream& input,T*& d)
    {bool data_exists;Read_Primitive(input,data_exists);
    if(data_exists){d=new T();Read_Binary<RW>(input,*d);}else d=0;} // potential memory leak if d pointed elsewhere

    static void Write(std::ostream& output,T* const& d)
    {Write_Primitive(output,d!=0); // write a bool tag indicating whether pointer's data follows
    if(d) Write_Binary<RW>(output,*d);}
};
//#####################################################################
// Read_Write for std::string's
//#####################################################################
template<class RW>
class Read_Write<std::string,RW>
{
public:
    static void Read(std::istream& input,std::string& d)
    {int n;Read_Primitive(input,n);
    char* buffer=new char[n];input.read(buffer,n);d.assign(buffer,buffer+n);delete[] buffer;}

    static void Write(std::ostream& output,const std::string& d)
    {int n=int(d.size());Write_Primitive(output,n);
    const char* s=d.c_str();output.write(s,n);}
};
//#####################################################################
}
#endif

#endif
