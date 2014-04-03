//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_VIEW
//#####################################################################
#ifndef __ARRAY_VIEW__
#define __ARRAY_VIEW__

#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<ARRAY_VIEW<T,ID> > {static const bool value=true;};
template<class T,class ID> struct IS_ARRAY_VIEW<ARRAY_VIEW<T,ID> > {static const bool value=true;};
template<class T,class ID,class T_NEW> struct REBIND<ARRAY_VIEW<T,ID>,T_NEW>{typedef ARRAY_VIEW<typename IF<IS_CONST<T>::value,const T_NEW,T_NEW>::TYPE,ID> TYPE;};

template<class T,class ID> struct IS_CONST<ARRAY_VIEW<const T,ID> > {static const bool value=true;}; // ARRAY_VIEW<const T,ID> is equivalent to const ARRAY_VIEW<const T,ID>

template<class T_ARRAY> struct CANONICALIZE_CONST_ARRAY<T_ARRAY,typename ENABLE_IF<IS_BASE_OF<ARRAY_VIEW<typename T_ARRAY::ELEMENT,typename T_ARRAY::INDEX>,T_ARRAY>::value>::TYPE>
    :public CANONICALIZE_CONST_ARRAY<ARRAY_VIEW<typename ADD_CONST<typename T_ARRAY::ELEMENT>::TYPE,typename T_ARRAY::INDEX> >{};

template<class T,class ID>
class ARRAY_VIEW:public ARRAY_BASE<typename REMOVE_CONST<T>::TYPE,ARRAY_VIEW<T,ID>,ID>
{
    struct UNUSABLE{};
    template<class S> struct COPY_CONST:public IF<IS_CONST<T>::value,typename ADD_CONST<S>::TYPE,S>{};
    typedef ARRAY_BASE<typename REMOVE_CONST<T>::TYPE,ARRAY_VIEW<T,ID>,ID> BASE;
public:
    typedef typename REMOVE_CONST<T>::TYPE ELEMENT;typedef ID INDEX;
    typedef T& RESULT_TYPE;

    // m and base_pointer inherit constness of T
    typename COPY_CONST<ID>::TYPE m;
    friend class ARRAY_VIEW<typename IF<IS_CONST<T>::value,ELEMENT,const ELEMENT>::TYPE,ID>;
    typename COPY_CONST<T*>::TYPE base_pointer;

public:
    using BASE::Same_Array;

    ARRAY_VIEW(const ID m,T* raw_data)
        :m(m),base_pointer(raw_data)
    {}

    ARRAY_VIEW(const ARRAY_VIEW<T,ID>& array)
        :m(array.m),base_pointer(array.base_pointer)
    {}

    ARRAY_VIEW(const ARRAY_VIEW<typename IF<IS_CONST<T>::value,typename REMOVE_CONST<T>::TYPE,UNUSABLE>::TYPE,ID>& array)
        :m(array.m),base_pointer(array.base_pointer)
    {}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY& array,typename ENABLE_IF<IS_SAME<ELEMENT,typename T_ARRAY::ELEMENT>::value && !IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
        :m(array.Size()),base_pointer(array.Get_Array_Pointer())
    {}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY array,typename ENABLE_IF<IS_SAME<ELEMENT,typename T_ARRAY::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
        :m(array.Size()),base_pointer(array.Get_Array_Pointer())
    {}

    ID Size() const
    {return m;}

    T& operator()(const ID i)
    {assert(ID(1)<=i && i<=m);return base_pointer[Value(i)-1];}

    const T& operator()(const ID i) const
    {assert(ID(1)<=i && i<=m);return base_pointer[Value(i)-1];}

    bool Valid_Index(const ID i) const
    {return ID(1)<=i && i<=m;}

    ARRAY_VIEW& operator=(const ARRAY_VIEW& source)
    {return BASE::operator=(source);}

    template<class T_ARRAY2>
    ARRAY_VIEW& operator=(const T_ARRAY2& source)
    {return BASE::operator=(source);}

    T* Get_Array_Pointer()
    {return base_pointer;}

    const T* Get_Array_Pointer() const
    {return base_pointer;}

    ARRAY_VIEW<typename REMOVE_CONST<T>::TYPE>& Const_Cast() const // return reference to allow Exchange
    {return reinterpret_cast<ARRAY_VIEW<typename REMOVE_CONST<T>::TYPE>&>(const_cast<ARRAY_VIEW&>(*this));}

    static bool Same_Array(const ARRAY_VIEW& array1,const ARRAY_VIEW& array2)
    {return array1.Get_Array_Pointer()==array2.Get_Array_Pointer();}

private:
    template<class T2>
    static void exchange(T2& a,T2& b)
    {T2 tmp=a;a=b;b=tmp;}

public:
    void Exchange(ARRAY_VIEW& other)
    {STATIC_ASSERT(!IS_CONST<T>::value); // make ARRAY_VIEW<const T> equivalent to const ARRAY_VIEW<const T>
    exchange(m,other.m);exchange(base_pointer,other.base_pointer);}

//#####################################################################
};
}
#endif
