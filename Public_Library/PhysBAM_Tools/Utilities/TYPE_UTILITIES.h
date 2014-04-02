//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TYPE_UTILITIES
//#####################################################################
#ifndef __TYPE_UTILITIES__
#define __TYPE_UTILITIES__

#include <cstdio>
#include <cstdlib>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>

namespace PhysBAM{
#define PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(trait_name,type_name,truth_value) \
template<> struct trait_name<type_name> {static const bool value=truth_value;}; \
template<> struct trait_name<const type_name> {static const bool value=truth_value;}; \
template<> struct trait_name<volatile type_name> {static const bool value=truth_value;}; \
template<> struct trait_name<const volatile type_name> {static const bool value=truth_value;};

template<bool b,class T=void> struct ENABLE_IF{};
template<class T> struct ENABLE_IF<true,T>{typedef T TYPE;};

template<bool b,class T=void> struct DISABLE_IF{};
template<class T> struct DISABLE_IF<false,T>{typedef T TYPE;};

template<bool b,class T1,class T2> struct IF{typedef T1 TYPE;};
template<class T1,class T2> struct IF<false,T1,T2>{typedef T2 TYPE;};

template<bool b> struct NOT {static const bool value=false;};
template<> struct NOT<false> {static const bool value=true;};

template<bool b1,bool b2,bool b3=true,bool b4=true> struct AND {static const bool value=b1 && b2 && b3 && b4;};
template<bool b1,bool b2,bool b3=false,bool b4=false> struct OR {static const bool value=b1 || b2 || b3 || b4;};

template<int i,int j> struct INTS_EQUAL {static const bool value=false;};
template<int i> struct INTS_EQUAL<i,i> {static const bool value=true;};

template<class T1,class T2> struct IS_SAME {static const bool value=false;};
template<class T> struct IS_SAME<T,T> {static const bool value=true;};

template<class T1,class T2> struct ASSERT_SAME_HELPER;
template<class T> struct ASSERT_SAME_HELPER<T,T>{};

#define STATIC_ASSERT_SAME(T1,T2) \
   typedef static_assert_test< \
      sizeof(::PhysBAM::ASSERT_SAME_HELPER<T1,T2>)> \
         PHYSBAM_JOIN(boost_static_assert_typedef_, __LINE__)

template<class T> struct IS_CONST {static const bool value=false;};
template<class T> struct IS_CONST<const T> {static const bool value=true;};

template<class T> struct REMOVE_CONST{typedef T TYPE;};
template<class T> struct REMOVE_CONST<const T>{typedef T TYPE;};

template<class T> struct ADD_CONST{typedef const T TYPE;};
template<class T> struct ADD_CONST<const T>{typedef T TYPE;};

template<class T> struct IS_POINTER {static const bool value=false;};
template<class T> struct IS_POINTER<T*> {static const bool value=true;};

template<class T> struct REMOVE_POINTER;
template<class T> struct REMOVE_POINTER<T*>{typedef T TYPE;};

template<class T> struct IS_MEMBER_POINTER {static const bool value=false;};
template<class T,class TF> struct IS_MEMBER_POINTER<TF T::*> {static const bool value=true;};

template<class T> struct IS_REFERENCE {static const bool value=false;};
template<class T> struct IS_REFERENCE<T&> {static const bool value=true;};

template<class T> struct ADD_REFERENCE{typedef T& TYPE;};
template<class T> struct ADD_REFERENCE<T&>{typedef T& TYPE;};

template<class T> struct REMOVE_REFERENCE{typedef T TYPE;};
template<class T> struct REMOVE_REFERENCE<T&>{typedef T TYPE;};

template<class T1,class T2=void,class T3=void,class T4=void> struct FIRST{typedef T1 TYPE;};

template<class T> struct IS_FLOAT_OR_DOUBLE:public OR<IS_SAME<T,float>::value,IS_SAME<T,double>::value>{};

template<class T_ARRAY,class ENABLER=void> struct IS_ARRAY {static const bool value=false;};
template<class T_ARRAY> struct IS_ARRAY<const T_ARRAY>:public IS_ARRAY<T_ARRAY>{};

template<class T_ARRAY,class ENABLER=void> struct IS_ARRAY_VIEW {static const bool value=false;};
template<class T_ARRAY> struct IS_ARRAY_VIEW<const T_ARRAY>:public IS_ARRAY_VIEW<T_ARRAY>{};

typedef char YES_TYPE;
struct NO_TYPE {char padding[8];};

template <class U> YES_TYPE Is_Class_Tester(void(U::*)(void));
template <class U> NO_TYPE Is_Class_Tester(...);

template <class T>
struct IS_CLASS
{
    static const bool union_value=false; // TODO: use __is_union(T) for gcc 4.3+
    static const bool value=AND<sizeof(Is_Class_Tester<T>(0)) == sizeof(YES_TYPE),NOT<union_value>::value>::value;
};

template<class T> struct IS_VOID {static const bool value=false;};
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_VOID,void,true)

template<class T> struct IS_INTEGRAL {static const bool value=false;};
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,unsigned char,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,unsigned short,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,unsigned int,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,unsigned long,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,unsigned long long,true)

PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,signed char,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,signed short,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,signed int,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,signed long,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,signed long long,true)

PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,bool,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,char,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_INTEGRAL,wchar_t,true)

template<class T> struct IS_FLOATING_POINT {static const bool value=false;};
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_FLOATING_POINT,float,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_FLOATING_POINT,double,true)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_FLOATING_POINT,long double,true)

template<class T> struct IS_SCALAR {static const bool value=OR<IS_INTEGRAL<T>::value,IS_FLOATING_POINT<T>::value>::value;};
template<class T> struct IS_FUNDAMENTAL {static const bool value=OR<IS_SCALAR<T>::value,IS_VOID<T>::value>::value;};

template <class T> struct REMOVE_CV_POINTER{};
template <class T> struct REMOVE_CV_POINTER<T*>{typedef T TYPE;};
template <class T> struct REMOVE_CV_POINTER<const T*>{typedef T TYPE;};
template <class T> struct REMOVE_CV_POINTER<volatile T*>{typedef T TYPE;};
template <class T> struct REMOVE_CV_POINTER<const volatile T*>{typedef T TYPE;};

template <class T> struct REMOVE_CV{typedef typename REMOVE_CV_POINTER<T*>::TYPE TYPE;};
template <class T> struct REMOVE_CV<T&>{typedef T& TYPE;};
template <class T,size_t N> struct REMOVE_CV<T const[N]>{typedef T TYPE[N];};
template <class T,size_t N> struct REMOVE_CV<T volatile[N]>{typedef T TYPE[N];};
template <class T,size_t N> struct REMOVE_CV<T const volatile[N]>{typedef T TYPE[N];};

template <class T> struct EMPTY_HELPER_BASE_T:public T{int i[256];};
struct EMPTY_HELPER_NOBASE{int i[256];};

template <class T,bool is_class=false> struct EMPTY_HELPER{static const bool value=false;};
template <class T> struct EMPTY_HELPER<T,true>{static const bool value=sizeof(EMPTY_HELPER_BASE_T<T>)==sizeof(EMPTY_HELPER_NOBASE);};

template <class T>
struct IS_EMPTY_VALUE
{
    typedef typename REMOVE_CV<T>::TYPE CV_TYPE;
    static const bool value=EMPTY_HELPER<CV_TYPE,IS_CLASS<CV_TYPE>::value>::value;
    // TODO: boost also ors this with __is_empty() for gcc 4.3+ to check for empty struct or union.
};
template <class T> struct IS_EMPTY {static const bool value=IS_EMPTY_VALUE<T>::value;};

// IS_BASE_OF
// from http://groups.google.com/group/comp.lang.c++.moderated/msg/dd6c4e4d5160bd83?pli=1
template <class B,class D>
struct IS_BASE_OF_VALUE
{
private:
    template<class T> static YES_TYPE check(D const volatile &, T);
    static NO_TYPE check(B const volatile &, int);

    struct C
    {
        operator B const volatile &() const;
        operator D const volatile &();
    };

    static C Get_C();
public:
    static const bool value=sizeof(check(Get_C(),0))==sizeof(YES_TYPE);
};
template <class B>
struct IS_BASE_OF_VALUE<B,B> {static const bool value=true;};
template <class B,class D> struct IS_BASE_OF {static const bool value=IS_BASE_OF_VALUE<B,D>::value;};

template <class T1, class T2>
struct IS_CONVERTIBLE_BASIC_IMPL
{
    static NO_TYPE _m_check(...);
    static YES_TYPE _m_check(T2);
    static T1 _m_from;

    static const bool value = sizeof(_m_check(_m_from))==sizeof(YES_TYPE);
};
template <class T1, class T2>
struct IS_CONVERTIBLE_IMPL
{
    typedef typename ADD_REFERENCE<T1>::TYPE ref_type;
    static const bool value=AND<IS_CONVERTIBLE_BASIC_IMPL<ref_type,T2>::value,NOT<IS_ARRAY<T2>::value>::value>::value;
};
struct TRUE_TYPE {static const bool value=true;};
struct FALSE_TYPE {static const bool value=false;};
template <bool trivial1, bool trivial2>
struct IS_CONVERTIBLE_IMPL_SELECT
{
   template <class T1, class T2>
   struct rebind
   {
      typedef IS_CONVERTIBLE_IMPL<T1,T2> type;
   };
};
template <>
struct IS_CONVERTIBLE_IMPL_SELECT<true, true>
{
   template <class T1, class T2>
   struct rebind
   {
      typedef TRUE_TYPE type;
   };
};
template<class T1,class T2> struct IS_CONVERTIBLE_DISPATCH_BASE
{
    typedef IS_CONVERTIBLE_IMPL_SELECT<IS_SCALAR<T1>::value,IS_SCALAR<T2>::value> selector;
    typedef typename selector::template rebind<T1,T2> isc_binder;
    typedef typename isc_binder::type type;
};
template<class T1,class T2> struct IS_CONVERTIBLE_DISPATCH:public IS_CONVERTIBLE_DISPATCH_BASE<T1,T2>::type{};
#define PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2_PART1(trait,spec1,spec2,base_value) \
    template<> struct trait<spec1,spec2> {static const bool value=base_value;}; \
    template<> struct trait<spec1,const spec2> {static const bool value=base_value;}; \
    template<> struct trait<spec1,volatile spec2> {static const bool value=base_value;}; \
    template<> struct trait<spec1,const volatile spec2> {static const bool value=base_value;};
#define PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2(trait,spec1,spec2,value) \
    PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2_PART1(trait,spec1,spec2,value) \
    PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2_PART1(trait,const spec1,spec2,value) \
    PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2_PART1(trait,volatile spec1,spec2,value) \
    PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2_PART1(trait,const volatile spec1,spec2,value)
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2(IS_CONVERTIBLE_IMPL,void,void,true)
#undef PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2
#undef PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS2_PART1
template<class T1,class T2> struct IS_CONVERTIBLE {static const bool value=IS_CONVERTIBLE_DISPATCH<T1,T2>::value;};

template<class T> struct IS_FUNCTION {static const bool value=NOT<IS_CONVERTIBLE<T*, const volatile void*>::value>::value;};
struct INT_CONVERTIBLE {INT_CONVERTIBLE(int);};
template <bool is_typename_arithmetic_or_reference = true> struct IS_ENUM_HELPER {template <class T> struct type {static const bool value=false;};};
template <> struct IS_ENUM_HELPER<false> {template <class T> struct type {static const bool value=IS_CONVERTIBLE<typename ADD_REFERENCE<T>::TYPE,INT_CONVERTIBLE>::value;};};
template <class T> struct IS_ENUM_IMPL
{
    static const bool selector=OR<OR<IS_SCALAR<T>::value,IS_REFERENCE<T>::value,IS_FUNCTION<T>::value,IS_CLASS<T>::value>::value,IS_ARRAY<T>::value>::value;
    typedef IS_ENUM_HELPER<selector> se_t;
    typedef typename se_t::template type<T> helper;
    static const bool value = helper::value;
};
PHYSBAM_TYPE_TRAIT_DECLARE_CV_QUALIFIERS(IS_ENUM_IMPL,void,false)
template<class T> struct IS_ENUM {static const bool value=IS_ENUM_IMPL<T>::value;};

template<class T,T v> struct STATIC_CONST {static const T value=v;};

template<class T> struct HAS_CHEAP_COPY:public OR<IS_FUNDAMENTAL<T>::value,IS_ENUM<T>::value,IS_ARRAY_VIEW<T>::value>{};

template<class T> struct REMOVE_ALL_EXTENTS {typedef T type;};
template<class T,int size> struct REMOVE_ALL_EXTENTS<T[size]> {typedef typename REMOVE_ALL_EXTENTS<T>::type type;};
template<class T> struct REMOVE_ALL_EXTENTS<T[]> {typedef typename REMOVE_ALL_EXTENTS<T>::type type;};

template<class T> struct IS_SCALAR_TYPE {static const bool value=IS_SCALAR<T>::value||IS_ENUM<T>::value||IS_POINTER<T>::value||IS_MEMBER_POINTER<T>::value;};

template<class T> struct IS_POD {static const bool value=IS_VOID<T>::value||IS_SCALAR_TYPE<typename REMOVE_ALL_EXTENTS<T>::type>::value;};

template<class T> struct HAS_TRIVIAL_DESTRUCTOR {static const bool value=IS_POD<T>::value;};

template<class T,class RW,class ENABLER=void> struct IS_BINARY_IO_SAFE;

template<class T,class SCALAR,class ENABLER=void> struct REPLACE_FLOATING_POINT{};
template<class T,class SCALAR> struct REPLACE_FLOATING_POINT<T,SCALAR,typename ENABLE_IF<AND<IS_FLOAT_OR_DOUBLE<T>::value,IS_FLOAT_OR_DOUBLE<SCALAR>::value>::value>::TYPE>{typedef SCALAR TYPE;};
template<class T,class SCALAR> struct REPLACE_FLOATING_POINT<T,SCALAR,typename ENABLE_IF<AND<AND<NOT<IS_FLOAT_OR_DOUBLE<T>::value>::value,IS_FUNDAMENTAL<T>::value>::value,IS_FLOAT_OR_DOUBLE<SCALAR>::value>::value>::TYPE>{typedef T TYPE;};
template<class T,class SCALAR> struct REPLACE_FLOATING_POINT<T,SCALAR,typename ENABLE_IF<IS_POINTER<T>::value>::TYPE> {typedef typename REPLACE_FLOATING_POINT<typename REMOVE_POINTER<T>::TYPE,SCALAR>::TYPE* TYPE;};
}
#endif
