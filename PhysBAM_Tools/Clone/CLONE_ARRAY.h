//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLONE_ARRAY 
//#####################################################################
#ifndef __CLONE_ARRAY__
#define __CLONE_ARRAY__

#include <PhysBAM_Tools/Clone/CLONEABLE.h>
namespace PhysBAM{

template<class T> class CLONE_ARRAY;

template<>
class CLONE_ARRAY<CLONEABLE_BASE>
{
private:
    const size_t sizeof_clone;
    const int count;
    char* data;
public:

    CLONE_ARRAY(const CLONEABLE_BASE& template_object,const int count); // constructs an array of default constructed type clones of the template object
    virtual ~CLONE_ARRAY();

    int Size() const
    {return count;}

    CLONEABLE_BASE& operator()(const int i)
    {return *(CLONEABLE_BASE*)(data+sizeof_clone*(i-1));}

    const CLONEABLE_BASE& operator()(const int i) const
    {return *(const CLONEABLE_BASE*)(data+sizeof_clone*(i-1));}
};

template<class T>
class CLONE_ARRAY:public CLONE_ARRAY<CLONEABLE_BASE>
{
    STATIC_ASSERT(IS_CLONEABLE<T>::value);
    typedef CLONE_ARRAY<CLONEABLE_BASE> BASE;
public:
    CLONE_ARRAY(const T& template_object,const int count)
        :CLONE_ARRAY<CLONEABLE_BASE>(template_object,count)
    {}

    T& operator()(const int i)
    {return static_cast<T&>(BASE::operator()(i));}

    const T& operator()(const int i) const
    {return static_cast<const T&>(BASE::operator()(i));}
};

//#####################################################################
}
#endif
