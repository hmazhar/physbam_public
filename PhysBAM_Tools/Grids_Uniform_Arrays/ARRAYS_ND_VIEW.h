//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_VIEW
//#####################################################################
#ifndef __ARRAYS_ND_VIEW__
#define __ARRAYS_ND_VIEW__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_BASE.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#endif
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,int d>
class ARRAY_VIEW<T,VECTOR<int,d> >:public ARRAYS_ND_BASE<VECTOR<typename REMOVE_CONST<T>::TYPE,d> >
{
    typedef VECTOR<int,d> TV_INT;
    struct UNUSABLE{};
    template<class S> struct COPY_CONST:public IF<IS_CONST<T>::value,typename ADD_CONST<S>::TYPE,S>{};
    typedef ARRAYS_ND_BASE<VECTOR<typename REMOVE_CONST<T>::TYPE,d> > BASE;
public:
    typedef typename REMOVE_CONST<T>::TYPE ELEMENT;typedef TV_INT INDEX;
    typedef T& RESULT_TYPE;

    using BASE::domain;using BASE::counts;using BASE::array;
private:
    friend class ARRAY_VIEW<typename IF<IS_CONST<T>::value,ELEMENT,const ELEMENT>::TYPE,TV_INT>;
    using BASE::base_pointer;

    using BASE::Calculate_Acceleration_Constants;

    void Initialize(typename REMOVE_CONST<T>::TYPE* raw_data)
    {int size=counts.Product();ARRAY_VIEW<typename REMOVE_CONST<T>::TYPE> new_array(size,raw_data);new_array.Exchange(array);Calculate_Acceleration_Constants();}

public:

    ARRAY_VIEW(const RANGE<TV_INT>& domain=RANGE<TV_INT>::Empty_Box(),typename REMOVE_CONST<T>::TYPE* raw_data=0)
        :BASE(domain)
    {Initialize(raw_data);}

    ARRAY_VIEW(const ARRAY_VIEW<typename REMOVE_CONST<T>::TYPE,TV_INT>& array)
        :BASE(array.domain)
    {Initialize(array.base_pointer);}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY& array,typename ENABLE_IF<IS_SAME<ELEMENT,typename T_ARRAY::ELEMENT>::value && !IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
        :BASE(array.Domain_Indices())
    {Initialize(array.Get_Array_Pointer());}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY array,typename ENABLE_IF<IS_SAME<ELEMENT,typename T_ARRAY::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
        :BASE(array.Domain_Indices())
    {Initialize(array.Get_Array_Pointer());}

    template<class T_ARRAY2>
    ARRAY_VIEW& operator=(const T_ARRAY2& source)
    {assert(Equal_Dimensions(*this,source));array=source.array;return *this;}

    T* Get_Array_Pointer()
    {return base_pointer;}

    const T* Get_Array_Pointer() const
    {return base_pointer;}

    void Exchange(ARRAY_VIEW& other)
    {STATIC_ASSERT(!IS_CONST<T>::value);this->Exchange(other);} // make ARRAY_VIEW<const T> equivalent to const ARRAY_VIEW<const T>

    static void Exchange_Arrays(ARRAY_VIEW& array1,ARRAY_VIEW& array2)
    {STATIC_ASSERT(!IS_CONST<T>::value); // make ARRAY_VIEW<const T> equivalent to const ARRAY_VIEW<const T>
    array1.array.Exchange(array2.array);
    exchange(array1.domain,array2.domain);exchange(array1.counts,array2.counts);array1.Calculate_Acceleration_Constants();array2.Calculate_Acceleration_Constants();}

    ARRAY_VIEW<typename REMOVE_CONST<T>::TYPE>& Const_Cast() const // return reference to allow Exchange
    {return reinterpret_cast<ARRAY_VIEW<typename REMOVE_CONST<T>::TYPE>&>(const_cast<ARRAY_VIEW&>(*this));}

    static bool Same_Array(const ARRAY_VIEW& array1,const ARRAY_VIEW& array2)
    {return array1.Get_Array_Pointer()==array2.Get_Array_Pointer();}
//#####################################################################
};
}
#endif
