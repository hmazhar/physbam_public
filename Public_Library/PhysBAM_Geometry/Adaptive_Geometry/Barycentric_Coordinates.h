//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GEOMETRIC_PREDICATES_DETAIL_BARYCENTRIC_COORDINATES__
#define __GEOMETRIC_PREDICATES_DETAIL_BARYCENTRIC_COORDINATES__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OBJECT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_SIGNED_VOLUME.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersection_Coordinates.h>
#include <cassert>

namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
//#################################################################
// Barycentric_Coordinates
//#################################################################
// computes adaptive coordinates
template<class T_EXACT,class T,int D,int N,class T_ADAPTIVE>
void Barycentric_Coordinates(const VECTOR<VECTOR<T,D>,D+1>& simplex,const VECTOR<VECTOR<T,D>,N>& point_vector,VECTOR<VECTOR<T_ADAPTIVE,D+1>,N>& barycentric_coordinates_vector)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,D>,D+1>,VECTOR<T,D> > T_RESULT;
    typedef typename T_RESULT::ADAPTIVE_TYPE1 ADAPTIVE_TYPE1;
    typedef typename T_RESULT::ADAPTIVE_TYPE2 ADAPTIVE_TYPE2;
    ADAPTIVE_TYPE2 common_denominator(Adaptive_Signed_Volume<T_EXACT>(simplex));
    VECTOR<ADAPTIVE_TYPE1,D+1> barycentric_coordinates_numerator;
    for(int i=1;i<=N;++i){
        Intersection_Coordinates_Helper<T_EXACT>(simplex,point_vector[i],barycentric_coordinates_numerator);
        for(int j=1;j<=D+1;++j)
            barycentric_coordinates_vector[i][j]=barycentric_coordinates_numerator[j]/common_denominator;}
}
//#################################################################
// Barycentric_Coordinates
//#################################################################
// computes non-adaptive coordinates to an absolute tolerance
template<class T_EXACT,class T,int D,int N,class U,class V>
void Barycentric_Coordinates(const VECTOR<VECTOR<T,D>,D+1>& simplex,const VECTOR<VECTOR<T,D>,N>& point_vector,VECTOR<VECTOR<U,D+1>,N>& barycentric_coordinates_vector,V abs_tol)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<U>::value));
    assert(abs_tol>0);
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,D>,D+1>,VECTOR<T,D> > T_RESULT;
    typedef typename T_RESULT::TYPE ADAPTIVE;
    VECTOR<VECTOR<ADAPTIVE,D+1>,N> adaptive_coordinates_vector;
    Barycentric_Coordinates<T_EXACT>(simplex,point_vector,adaptive_coordinates_vector);
    for(int i=1;i<=N;++i)
        Convert_And_Ensure_Abs_Tol(adaptive_coordinates_vector[i],barycentric_coordinates_vector[i],abs_tol);
}
//#################################################################
}
using GEOMETRIC_PREDICATES_DETAIL::Barycentric_Coordinates;
}
#endif
