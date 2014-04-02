//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Find_Type
//#####################################################################
#ifndef __Find_Type__
#define __Find_Type__

#include <typeinfo>
namespace PhysBAM{

template<class T> struct FIND_TYPE_HELPER; // T must be a pointer or a reference

template<class T> struct FIND_TYPE_HELPER<T*>
{
    template<class T_ARRAY> static T* Find_Type(T_ARRAY& array,const int index)
    {int found=0;
    for(typename T_ARRAY::INDEX i(1);i<=array.Size();i++){
        T* element=dynamic_cast<T*>(array(i));
        if(element){found++;if(found==index) return element;}}
    return 0;}
};

template<class T> struct FIND_TYPE_HELPER<T&>
{
    template<class T_ARRAY> static T& Find_Type(const T_ARRAY& array,const int index)
    {T* element=FIND_TYPE_HELPER<T*>::Find_Type(array,index);
    if(!element) throw std::bad_cast();
    return *element;}
};

template<class T,class T_ARRAY> T Find_Type(T_ARRAY& array,const int index=1)
{
    return FIND_TYPE_HELPER<T>::Find_Type(array,index);
}

}
#endif
