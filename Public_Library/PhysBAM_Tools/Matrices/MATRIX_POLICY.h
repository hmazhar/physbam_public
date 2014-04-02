//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_POLICY
//#####################################################################
#ifndef __MATRIX_POLICY__
#define __MATRIX_POLICY__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV> struct MATRIX_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct MATRIX_POLICY<VECTOR<T,0> >
{
    typedef MATRIX<T,0> SYMMETRIC_MATRIX;
    typedef MATRIX<T,0> DIAGONAL_MATRIX;
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct MATRIX_POLICY<VECTOR<T,1> >
{
    typedef MATRIX<T,1> SYMMETRIC_MATRIX;
    typedef MATRIX<T,1> DIAGONAL_MATRIX;
    typedef PhysBAM::MATRIX<T,2> TRANSFORMATION_MATRIX;
    typedef PhysBAM::MATRIX<T,0,1> CROSS_PRODUCT_MATRIX;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct MATRIX_POLICY<VECTOR<T,2> >
{
    typedef PhysBAM::SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX;
    typedef PhysBAM::DIAGONAL_MATRIX<T,2> DIAGONAL_MATRIX;
    typedef PhysBAM::MATRIX<T,3> TRANSFORMATION_MATRIX;
    typedef PhysBAM::MATRIX<T,1,2> CROSS_PRODUCT_MATRIX;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct MATRIX_POLICY<VECTOR<T,3> >
{
    typedef PhysBAM::SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX;
    typedef PhysBAM::DIAGONAL_MATRIX<T,3> DIAGONAL_MATRIX;
    typedef PhysBAM::MATRIX<T,4> TRANSFORMATION_MATRIX;
    typedef PhysBAM::MATRIX<T,3> CROSS_PRODUCT_MATRIX;
};
//#####################################################################
}
#endif
