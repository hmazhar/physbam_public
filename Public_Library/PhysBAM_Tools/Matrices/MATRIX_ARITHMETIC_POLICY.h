//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_ARITHMETIC_POLICY
//#####################################################################
#ifndef __MATRIX_ARITHMETIC_POLICY__
#define __MATRIX_ARITHMETIC_POLICY__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,VECTOR<T,n> >{typedef VECTOR<T,m> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,VECTOR<T,d> >{typedef VECTOR_ND<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};
template<class T> struct PRODUCT<MATRIX_MXN<T>,VECTOR_ND<T> >{typedef VECTOR_ND<T> TYPE;};

template<class T,int m,int n,int k> struct PRODUCT<MATRIX<T,m,k>,MATRIX<T,k,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,DIAGONAL_MATRIX<T,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,SYMMETRIC_MATRIX<T,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,UPPER_TRIANGULAR_MATRIX<T,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int m,int n> struct PRODUCT<MATRIX<T,m,n>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<DIAGONAL_MATRIX<T,m>,MATRIX<T,m,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,DIAGONAL_MATRIX<T,d> >{typedef DIAGONAL_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef UPPER_TRIANGULAR_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<DIAGONAL_MATRIX<T,d>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<SYMMETRIC_MATRIX<T,m>,MATRIX<T,m,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,DIAGONAL_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<SYMMETRIC_MATRIX<T,d>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,m>,MATRIX<T,m,n> >{typedef MATRIX<T,m,n> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,DIAGONAL_MATRIX<T,d> >{typedef UPPER_TRIANGULAR_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef UPPER_TRIANGULAR_MATRIX<T,d> TYPE;};
template<class T,int d> struct PRODUCT<UPPER_TRIANGULAR_MATRIX<T,d>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T,int m,int n> struct PRODUCT<MATRIX_MXN<T>,MATRIX<T,m,n> >{typedef MATRIX_MXN<T> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,DIAGONAL_MATRIX<T,d> >{typedef MATRIX_MXN<T> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,SYMMETRIC_MATRIX<T,d> >{typedef MATRIX_MXN<T> TYPE;};
template<class T,int d> struct PRODUCT<MATRIX_MXN<T>,UPPER_TRIANGULAR_MATRIX<T,d> >{typedef MATRIX_MXN<T> TYPE;};
template<class T> struct PRODUCT<MATRIX_MXN<T>,MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class TV> struct PRODUCT<FRAME<TV>,TV>{typedef TV TYPE;};

template<class TM> struct TRANSPOSE;
template<class T,int m,int n> struct TRANSPOSE<MATRIX<T,m,n> >{typedef MATRIX<T,n,m> TYPE;};
template<class T,int d> struct TRANSPOSE<SYMMETRIC_MATRIX<T,d> >{typedef SYMMETRIC_MATRIX<T,d> TYPE;};
template<class T,int d> struct TRANSPOSE<DIAGONAL_MATRIX<T,d> >{typedef DIAGONAL_MATRIX<T,d> TYPE;};
template<class T,int d> struct TRANSPOSE<UPPER_TRIANGULAR_MATRIX<T,d> >{typedef MATRIX<T,d> TYPE;};
template<class T> struct TRANSPOSE<MATRIX_MXN<T> >{typedef MATRIX_MXN<T> TYPE;};

template<class T1,class T2> struct PRODUCT_TRANSPOSE{typedef typename PRODUCT<T1,typename TRANSPOSE<T2>::TYPE>::TYPE TYPE;};
template<class T1,class T2> struct TRANSPOSE_PRODUCT{typedef typename PRODUCT<typename TRANSPOSE<T1>::TYPE,T2>::TYPE TYPE;};

}
#endif
