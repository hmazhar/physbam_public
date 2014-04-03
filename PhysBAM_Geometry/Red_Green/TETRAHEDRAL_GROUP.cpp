//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAEHDRAL_GROUP
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Red_Green/TETRAHEDRAL_GROUP.h>
using namespace PhysBAM;
//#####################################################################
// source: http://www.math.sunysb.edu/~tony/bintet
template<class T> const ROTATION<VECTOR<T,3> > TETRAHEDRAL_GROUP<T>::
quaternions[12]={ROTATION<VECTOR<T,3> >::From_Components(1,0,0,0),ROTATION<VECTOR<T,3> >::From_Components(0,1,0,0),ROTATION<VECTOR<T,3> >::From_Components(0,0,1,0),ROTATION<VECTOR<T,3> >::From_Components(0,0,0,1),ROTATION<VECTOR<T,3> >::From_Components(1,-1,-1,-1),ROTATION<VECTOR<T,3> >::From_Components(-1,-1,-1,-1),
    ROTATION<VECTOR<T,3> >::From_Components(1,-1,1,1),ROTATION<VECTOR<T,3> >::From_Components(-1,-1,1,1),ROTATION<VECTOR<T,3> >::From_Components(1,1,-1,1),ROTATION<VECTOR<T,3> >::From_Components(-1,1,-1,1),ROTATION<VECTOR<T,3> >::From_Components(1,1,1,-1),ROTATION<VECTOR<T,3> >::From_Components(-1,1,1,-1)};
template<class T> const ROTATION<VECTOR<int,3> > TETRAHEDRAL_GROUP<T>::
quaternions_times_2[12]={ROTATION<VECTOR<int,3> >::From_Components(2,0,0,0),ROTATION<VECTOR<int,3> >::From_Components(0,2,0,0),ROTATION<VECTOR<int,3> >::From_Components(0,0,2,0),ROTATION<VECTOR<int,3> >::From_Components(0,0,0,2),ROTATION<VECTOR<int,3> >::From_Components(1,-1,-1,-1),ROTATION<VECTOR<int,3> >::From_Components(-1,-1,-1,-1),
    ROTATION<VECTOR<int,3> >::From_Components(1,-1,1,1),ROTATION<VECTOR<int,3> >::From_Components(-1,-1,1,1),ROTATION<VECTOR<int,3> >::From_Components(1,1,-1,1),ROTATION<VECTOR<int,3> >::From_Components(-1,1,-1,1),ROTATION<VECTOR<int,3> >::From_Components(1,1,1,-1),ROTATION<VECTOR<int,3> >::From_Components(-1,1,1,-1)};
template<class T> const int TETRAHEDRAL_GROUP<T>::multiplication_table[12][12]={
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11},
    { 1, 0, 3, 2,10, 9, 8,11, 6, 5, 4, 7},
    { 2, 3, 0, 1, 6,11, 4, 9,10, 7, 8, 5},
    { 3, 2, 1, 0, 8, 7,10, 5, 4,11, 6, 9},
    { 4, 8,10, 6, 5, 0, 9, 2,11, 3, 7, 1},
    { 5,11, 7, 9, 0, 4, 3,10, 1, 6, 2, 8},
    { 6,10, 8, 4,11, 2, 7, 0, 5, 1, 9, 3},
    { 7, 9, 5,11, 3, 8, 0, 6, 2,10, 1, 4},
    { 8, 4, 6,10, 7, 3,11, 1, 9, 0, 5, 2},
    { 9, 7,11, 5, 1,10, 2, 4, 0, 8, 3, 6},
    {10, 6, 4, 8, 9, 1, 5, 3, 7, 2,11, 0},
    {11, 5, 9, 7, 2, 6, 1, 8, 3, 4, 0,10}};
template<class T> const int TETRAHEDRAL_GROUP<T>::inversion_table[12]={0,1,2,3,5,4,7,6,9,8,11,10};
template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::e=0;template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::i=1;
template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::j=2;template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::k=3;
template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::a=4;template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::a2=5;
template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::b=6;template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::b2=7;
template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::c=8;template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::c2=9;
template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::d=10;template<class T> const TETRAHEDRAL_GROUP<T> TETRAHEDRAL_GROUP<T>::d2=11;
//#####################################################################
template class TETRAHEDRAL_GROUP<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TETRAHEDRAL_GROUP<double>;
#endif
#endif
