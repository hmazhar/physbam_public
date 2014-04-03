//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARITHMETIC_POLICY
//#####################################################################
#ifndef __ARITHMETIC_POLICY__
#define __ARITHMETIC_POLICY__

#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T1,class T2,class ENABLER=void> struct CAN_ASSIGN;
template<class T1,class T2,class ENABLER=void> struct SUM;
template<class T1,class T2,class ENABLER=void> struct DIFFERENCE;
template<class T1,class T2,class ENABLER=void> struct PRODUCT;
template<class T1,class T2,class ENABLER=void> struct QUOTIENT;
template<class T,class ENABLER=void> struct NEGATION;

// Builtin
template<class T> struct CAN_ASSIGN<T,T>{static const bool value=true;};
template<> struct CAN_ASSIGN<double,float>{static const bool value=false;};
template<> struct CAN_ASSIGN<double,int>{static const bool value=false;};
template<> struct CAN_ASSIGN<float,double>{static const bool value=false;};
template<> struct CAN_ASSIGN<float,int>{static const bool value=false;};
template<> struct CAN_ASSIGN<int,double>{static const bool value=false;};
template<> struct CAN_ASSIGN<int,float>{static const bool value=false;};

template<> struct SUM<double,double>{typedef double TYPE;};
template<> struct SUM<double,float>{typedef double TYPE;};
template<> struct SUM<double,int>{typedef double TYPE;};
template<> struct SUM<float,double>{typedef double TYPE;};
template<> struct SUM<float,float>{typedef float TYPE;};
template<> struct SUM<float,int>{typedef float TYPE;};
template<> struct SUM<int,double>{typedef double TYPE;};
template<> struct SUM<int,float>{typedef float TYPE;};
template<> struct SUM<int,int>{typedef int TYPE;};

template<> struct DIFFERENCE<double,double>{typedef double TYPE;};
template<> struct DIFFERENCE<double,float>{typedef double TYPE;};
template<> struct DIFFERENCE<double,int>{typedef double TYPE;};
template<> struct DIFFERENCE<float,double>{typedef double TYPE;};
template<> struct DIFFERENCE<float,float>{typedef float TYPE;};
template<> struct DIFFERENCE<float,int>{typedef float TYPE;};
template<> struct DIFFERENCE<int,double>{typedef double TYPE;};
template<> struct DIFFERENCE<int,float>{typedef float TYPE;};
template<> struct DIFFERENCE<int,int>{typedef int TYPE;};

template<> struct PRODUCT<double,double>{typedef double TYPE;};
template<> struct PRODUCT<double,float>{typedef double TYPE;};
template<> struct PRODUCT<double,int>{typedef double TYPE;};
template<> struct PRODUCT<float,double>{typedef double TYPE;};
template<> struct PRODUCT<float,float>{typedef float TYPE;};
template<> struct PRODUCT<float,int>{typedef float TYPE;};
template<> struct PRODUCT<int,double>{typedef double TYPE;};
template<> struct PRODUCT<int,float>{typedef float TYPE;};
template<> struct PRODUCT<int,int>{typedef int TYPE;};

template<> struct QUOTIENT<double,double>{typedef double TYPE;};
template<> struct QUOTIENT<double,float>{typedef double TYPE;};
template<> struct QUOTIENT<double,int>{typedef double TYPE;};
template<> struct QUOTIENT<float,double>{typedef double TYPE;};
template<> struct QUOTIENT<float,float>{typedef float TYPE;};
template<> struct QUOTIENT<float,int>{typedef float TYPE;};
template<> struct QUOTIENT<int,double>{typedef double TYPE;};
template<> struct QUOTIENT<int,float>{typedef float TYPE;};
template<> struct QUOTIENT<int,int>{typedef int TYPE;};

template<> struct NEGATION<double>{typedef double TYPE;};
template<> struct NEGATION<float>{typedef float TYPE;};
template<> struct NEGATION<int>{typedef int TYPE;};

//template<class T,class TV> struct SUM<T,TV,typename ENABLE_IF<AND<IS_SCALAR<T>::value,NOT<IS_SCALAR<TV>::value>::value>::value>::TYPE>{typedef TV TYPE;};
//template<class T,class TV> struct SUM<TV,T,typename ENABLE_IF<AND<IS_SCALAR<T>::value,NOT<IS_SCALAR<TV>::value>::value>::value>::TYPE>{typedef TV TYPE;};
//template<class T,class TV> struct PRODUCT<T,TV,typename ENABLE_IF<IS_SCALAR<T>::value>::TYPE>{typedef TV TYPE;};

}
#endif
