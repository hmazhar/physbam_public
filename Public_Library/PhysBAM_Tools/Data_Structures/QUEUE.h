//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUEUE
//#####################################################################
#ifndef __QUEUE__
#define __QUEUE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

template<class T>
class QUEUE
{
public:
    ARRAY<T> array;
    int front,back;

    explicit QUEUE(const int size)
         :array(size+1),front(1),back(1)
    {}

    void Enqueue(const T& element)
    {array(back)=element;if(++back>array.m) back=1;
    assert(back!=front);} // dies if you run out of room

    void Safe_Enqueue(const T& element) // resizes array on overflow
    {array(back)=element;if(++back>array.m) back=1;
    if(back==front){
        ARRAY<T> new_array(4*array.m/3+2);back=1;
        for(int index=front;index<=array.m;index++) new_array(back++)=array(index);
        for(int index=1;index<front;index++) new_array(back++)=array(index);
        front=1;ARRAY<T>::Exchange_Arrays(array,new_array);}}

    T Dequeue()
    {assert(!Empty());int index=front;if(++front>array.m) front=1;return array(index);}

    const T& Peek() const
    {return array(front);}

    const T& operator()(const int i)
    {assert(i<=Size());int index=front+i-1;if(index>array.m) index=index-array.m;return array(index);}

    int Size() const
    {if(back<front) return back+array.m-front;else return back-front;}

    bool Empty() const
    {return back==front;}

    bool Full() const
    {return Size()==array.m-1;}

    void Remove_All()
    {front=1;back=1;}

//#####################################################################
};
}
#endif
