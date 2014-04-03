//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STACK
//##################################################################### 
#ifndef __STACK__
#define __STACK__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

template<class T>
class STACK
{
public:
    ARRAY<T> array;

    STACK()
    {}

    void Clean_Memory()
    {array.Clean_Memory();}

    void Preallocate(const int max_size)
    {array.Preallocate(max_size);}

    void Increase_Size(const int size)
    {array.Preallocate(size+array.m);}

    void Compact()
    {array.Compact();}

    void Push(const T& element) PHYSBAM_ALWAYS_INLINE
    {array.Append(element);}

    T Pop()
    {return array.Pop();}
 
    const T& Peek() const
    {return array.Last();}

    bool Empty() const
    {return array.m==0;}

    void Remove_All()
    {array.Remove_All();}

    void Exchange(STACK& stack)
    {array.Exchange(stack.array);}

//#####################################################################
};  
}
#endif
