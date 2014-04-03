//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header MATRIX_FORWARD
//#####################################################################
#ifndef __MATRIX_FORWARD__
#define __MATRIX_FORWARD__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class T_MATRIX> class MATRIX_BASE;
template<class T,int m_input,int n_input=m_input> class MATRIX;
template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d> class SYMMETRIC_MATRIX;
template<class T,int d> class UPPER_TRIANGULAR_MATRIX;
template<class T> class MATRIX_MXN;
template<class T> class SYMMETRIC_MATRIX_NXN;
template<class T> class SPARSE_MATRIX_NXN;
template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T_MATRIX> class TRANSPOSE_MATRIX;

template<class T> struct IS_SCALAR_BLOCK;
template<class T> struct IS_SCALAR_VECTOR_SPACE;
template<class T,int m,int n> struct IS_SCALAR_BLOCK<MATRIX<T,m,n> > {static const bool value=IS_SCALAR_BLOCK<T>::value;};
template<class T,int m,int n> struct IS_SCALAR_VECTOR_SPACE<MATRIX<T,m,n> > {static const bool value=IS_SCALAR_VECTOR_SPACE<T>::value;};
template<class T,int m,int n,class RW> struct IS_BINARY_IO_SAFE<MATRIX<T,m,n>,RW> {static const bool value=AND<(m>0),(n>0),IS_BINARY_IO_SAFE<T,RW>::value>::value;};
template<class T,int m,int n,class SCALAR> struct REPLACE_FLOATING_POINT<MATRIX<T,m,n>,SCALAR> {typedef MATRIX<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,m,n> TYPE;};
template<class T,int d,class SCALAR> struct REPLACE_FLOATING_POINT<DIAGONAL_MATRIX<T,d>,SCALAR> {typedef DIAGONAL_MATRIX<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,d> TYPE;};

}
#endif
