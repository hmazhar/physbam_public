//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Functions Sort and Stable_Sort
//#####################################################################
#ifndef __SORT__
#define __SORT__

#include <PhysBAM_Tools/Math_Tools/min.h>
#include <algorithm>
#include <functional>
namespace PhysBAM{

//#####################################################################
// Comparison Function Objects
//#####################################################################
template<class T_ARRAY,class T_COMPARE=std::less<typename T_ARRAY::ELEMENT> >
struct INDIRECT_COMPARE
{
    const T_ARRAY& array;
    T_COMPARE comparison;

    INDIRECT_COMPARE(const T_ARRAY& array_input)
        :array(array_input)
    {}

    INDIRECT_COMPARE(const T_ARRAY& array_input,const T_COMPARE& comparison_input)
        :array(array_input),comparison(comparison_input)
    {}

    bool operator()(const int index1,const int index2) const
    {return comparison(array(index1),array(index2));}
};

template<class T_ARRAY>
inline INDIRECT_COMPARE<T_ARRAY> Indirect_Comparison(const T_ARRAY& array)
{return INDIRECT_COMPARE<T_ARRAY>(array);}

template<class T_ARRAY,class T_COMPARISON>
inline INDIRECT_COMPARE<T_ARRAY,T_COMPARISON> Indirect_Comparison(const T_ARRAY& array,const T_COMPARISON& comparison)
{return INDIRECT_COMPARE<T_ARRAY,T_COMPARISON>(array,comparison);}

template<class T,class T_FIELD>
struct FIELD_COMPARE
{
    T_FIELD T::*field;

    FIELD_COMPARE(T_FIELD T::*field_input)
        :field(field_input)
    {}

    bool operator()(const T& x1,const T& x2) const
    {return x1.*field<x2.*field;}
};

template<class T,class T_FIELD>
inline FIELD_COMPARE<T,T_FIELD> Field_Comparison(T_FIELD T::*field)
{return FIELD_COMPARE<T,T_FIELD>(field);}

struct LEXICOGRAPHIC_COMPARE
{
    template<class T_ARRAY> bool operator()(const T_ARRAY& a1,const T_ARRAY& a2) const
    {int m1=a1.Size(),m2=a2.Size(),m=min(m1,m2);
    for(int i=1;i<=m;i++) if(a1(i)!=a2(i)) return a1(i)<a2(i);
    return m1<m2;}
};

//#####################################################################
// Functions Sort and Stable_Sort
//#####################################################################
template<class T_ARRAY,class T_COMPARE>
inline void Sort(T_ARRAY& array,const T_COMPARE& comparison)
{
    std::sort(array.begin(),array.end(),comparison);
}

template<class T_ARRAY,class T_COMPARE>
inline void Stable_Sort(T_ARRAY& array,const T_COMPARE& comparison)
{
    std::stable_sort(array.begin(),array.end(),comparison);
}

template<class T_ARRAY>
inline void Sort(T_ARRAY& array)
{Sort(array,std::less<typename T_ARRAY::ELEMENT>());}

template<class T_ARRAY>
inline void Stable_Sort(T_ARRAY& array)
{Stable_Sort(array,std::less<typename T_ARRAY::ELEMENT>());}

//#####################################################################
}
#endif
