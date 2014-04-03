//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace ARRAYS_COMPUTATIONS
//#####################################################################
#ifndef __ARRAY_HEAPIFY__
#define __ARRAY_HEAPIFY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
namespace PhysBAM{

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T2,class ID>
    void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux,ID index,const ID heap_size) // largest on top, only sorts down from index (not up!)
    {for(;;){ID left(2*Value(index)),right(2*Value(index)+1),index_of_largest=index;
        if(left<=heap_size && a(left)>a(index_of_largest)) index_of_largest=left;
        if(right<=heap_size && a(right)>a(index_of_largest)) index_of_largest=right;
        if(index_of_largest!=index){exchange(a(index),a(index_of_largest));exchange(aux(index),aux(index_of_largest));index=index_of_largest;}else return;}}

    template<class T,class T2,class ID>
    void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux) // largest on top
    {for(ID i=a.m/2;i>=ID(1);i--) Heapify(a,aux,i,a.m);}

    template<class T,class T2,class ID>
    void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux,const ID max_index) // largest on top, only does from 1 to max_index
    {for(ID i(Value(max_index/2));i>=ID(1);i--) Heapify(a,aux,i,max_index);}

    template<class T_ARRAY>
    void Heapify(T_ARRAY& a,typename T_ARRAY::INDEX index,const typename T_ARRAY::INDEX heap_size) // largest on top, only sorts down from index (not up!)
    {typedef typename T_ARRAY::INDEX ID;
    for(;;){ID left(2*Value(index)),right(2*Value(index)+1),index_of_largest=index;
        if(left<=heap_size && a(left)>a(index_of_largest)) index_of_largest=left;
        if(right<=heap_size && a(right)>a(index_of_largest)) index_of_largest=right;
        if(index_of_largest!=index){exchange(a(index),a(index_of_largest));index=index_of_largest;}else return;}}

    template<class T_ARRAY>
    void Reverse_In_Place(T_ARRAY& input)
    {typedef typename T_ARRAY::INDEX ID;for(ID i(1);i<=ID(Value(input.m)/2);i++) exchange(input(i),input(input.m+1-i));}

    template<class T_ARRAY>
    void Heapify(T_ARRAY& a) // largest on top
    {typedef typename T_ARRAY::INDEX ID;for(ID i=a.Size()/2;i>=ID(1);i--) Heapify(a,i,a.Size());}

    template<class T,class ID>
    void Heapify(ARRAY<T,ID>& a,const ID max_index) // largest on top, only does from 1 to max_index
    {for(ID i=max_index/2;i>=ID(1);i--) Heapify(a,i,max_index);}

    template<class T,class ID>
    void Compact_Array_Using_Compaction_Array(ARRAY<T,ID>& array,const ARRAY<ID,ID>& compaction_array,ARRAY<T,ID>* temporary_array=0)
    {ID compaction_array_m=compaction_array.Size();
    bool temporary_array_defined=temporary_array!=0;if(!temporary_array_defined) temporary_array=new ARRAY<T,ID>(compaction_array_m,false);
    ARRAY<T,ID>::Put(array,*temporary_array);for(ID i(1);i<=compaction_array_m;i++) if(compaction_array(i)>0) array(compaction_array(i))=(*temporary_array)(i);
    if(!temporary_array_defined){delete temporary_array;temporary_array=0;}}

//#####################################################################
}
}
#endif
