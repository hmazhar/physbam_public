//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERSECTION_COORDINATES__
#define __INTERSECTION_COORDINATES__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OBJECT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_QUOTIENT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_SIGNED_VOLUME.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_SUM.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <cassert>

namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
//#################################################################
// class INTERSECTION_COORDINATES_ADAPTIVE_RESULT
//#################################################################
// contains typedefs EXACT_TYPE and TYPE to easily deduce a proper adaptive type to receive intersection coordinates depending on
//  the types of the input arguments passed to Intersection_Coordinates i.e. TYPE == ADAPTIVE_OBJECT<EXACT_TYPE>
// TODO: Perhaps move this into a namespace result_of, as in boost?

template<class T0,class T1,class T2,class T3=void>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT;

template<class T0,class T1,class T2,class T3=void>
struct INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT;

template<class EXACT_TYPE,class T,int N>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE,VECTOR<T,N>,VECTOR<VECTOR<T,N>,N+1> >
    :public INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE,VECTOR<VECTOR<T,N>,N+1>,VECTOR<T,N> >
{};

template<class EXACT_TYPE_T,class T,int N>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,N>,N+1>,VECTOR<T,N> >
{
    typedef typename INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,N>,N+1>,VECTOR<T,N> >::EXACT_TYPE EXACT_TYPE1;
    typedef ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<EXACT_TYPE_T,T>::ADAPTIVE,N> SIGNED_VOLUME_RESULT;
    typedef typename SIGNED_VOLUME_RESULT::EXACT_TYPE EXACT_TYPE2;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE1> ADAPTIVE_TYPE1;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE2> ADAPTIVE_TYPE2;
    typedef typename ADAPTIVE_QUOTIENT<ADAPTIVE_TYPE1,ADAPTIVE_TYPE2>::EXACT_TYPE EXACT_TYPE;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE> TYPE;
};

template<class EXACT_TYPE_T,class T>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,2>,2>,VECTOR<VECTOR<T,2>,2> >
{
    typedef ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<EXACT_TYPE_T,T>::ADAPTIVE,2> SIGNED_VOLUME_RESULT;
    typedef typename SIGNED_VOLUME_RESULT::EXACT_TYPE EXACT_TYPE1;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE1> ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_SUM<ADAPTIVE_TYPE1,ADAPTIVE_TYPE1>::EXACT_TYPE EXACT_TYPE2;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE2> ADAPTIVE_TYPE2;
    typedef typename ADAPTIVE_QUOTIENT<ADAPTIVE_TYPE1,ADAPTIVE_TYPE2>::EXACT_TYPE EXACT_TYPE;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE> TYPE;
};

template<class EXACT_TYPE,class T>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE,VECTOR<VECTOR<T,3>,2>,VECTOR<VECTOR<T,3>,3> >
    :public INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,2> >
{};

template<class EXACT_TYPE_T,class T>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,2> >
{
    typedef ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<EXACT_TYPE_T,T>::ADAPTIVE,3> SIGNED_VOLUME_RESULT;
    typedef typename SIGNED_VOLUME_RESULT::EXACT_TYPE EXACT_TYPE1;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE1> ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_SUM<ADAPTIVE_TYPE1,ADAPTIVE_TYPE1>::EXACT_TYPE EXACT_TYPE2;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE2> ADAPTIVE_TYPE2;
    typedef typename ADAPTIVE_QUOTIENT<ADAPTIVE_TYPE1,ADAPTIVE_TYPE2>::EXACT_TYPE EXACT_TYPE;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE> TYPE;
};

template<class EXACT_TYPE_T,class T>
struct INTERSECTION_COORDINATES_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> >
{
    typedef INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::TYPE2 ADAPTIVE_TYPE2;
    typedef typename ADAPTIVE_QUOTIENT<ADAPTIVE_TYPE1,ADAPTIVE_TYPE2>::EXACT_TYPE EXACT_TYPE;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE> TYPE;
};

template<class EXACT_TYPE_T,class T,int N>
struct INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,N>,N+1>,VECTOR<T,N> >
{
    typedef ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<EXACT_TYPE_T,T>::ADAPTIVE,N> SIGNED_VOLUME_RESULT;
    typedef typename SIGNED_VOLUME_RESULT::EXACT_TYPE EXACT_TYPE;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE> TYPE;
};

template<class EXACT_TYPE_T,class T>
struct INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<EXACT_TYPE_T,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> >
{
    typedef ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<EXACT_TYPE_T,T>::ADAPTIVE,3> SIGNED_VOLUME_RESULT;
    typedef typename SIGNED_VOLUME_RESULT::EXACT_TYPE EXACT_TYPE0;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE0> ADAPTIVE_TYPE0;
    typedef ADAPTIVE_DETERMINANT<ADAPTIVE_TYPE0,2> DETERMINANT_RESULT;
    typedef typename DETERMINANT_RESULT::EXACT_TYPE EXACT_TYPE1;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE1> TYPE1;
    typedef typename ADAPTIVE_SUM_TYPE<TYPE1,TYPE1,TYPE1>::TYPE::EXACT_TYPE EXACT_TYPE2;
    typedef ADAPTIVE_OBJECT<EXACT_TYPE2> TYPE2;
};
//#################################################################
// Intersection_Coordinates_Helper
//#################################################################
template<class T_EXACT,class T,class T_ADAPTIVE>
void Intersection_Coordinates_Helper(const VECTOR<VECTOR<T,2>,3>& triangle,const VECTOR<T,2>& point,VECTOR<T_ADAPTIVE,3>& triangle_coordinates_numerator)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    // triangle = [p1 p2 p3]
    // point = q
    // x = a1 * p1 + a2 * p2 + a3 * p3
    //   = q
    // 1 = a1 + a2 + a3
    // Da1 = [ q p2 p3 ]
    // Da2 = [ p1 q p3 ]
    // Da3 = [ p1 p2 q ]
    triangle_coordinates_numerator[1]=Adaptive_Signed_Volume<T_EXACT>(point,triangle[2],triangle[3]);
    triangle_coordinates_numerator[2]=Adaptive_Signed_Volume<T_EXACT>(triangle[1],point,triangle[3]);
    triangle_coordinates_numerator[3]=Adaptive_Signed_Volume<T_EXACT>(triangle[1],triangle[2],point);
    // Da  = [ p1 p2 p3 ]
    // a1  = Da1 / Da
    // a2  = Da2 / Da
    // a3  = Da3 / Da
}

template<class T_EXACT,class T,class T_ADAPTIVE>
void Intersection_Coordinates_Helper(const VECTOR<VECTOR<T,3>,4>& tetrahedron,const VECTOR<T,3>& point,VECTOR<T_ADAPTIVE,4>& tetrahedron_coordinates_numerator)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    // tetrahedron = [p1 p2 p3 p4]
    // point = q
    // x = a1 * p1 + a2 * p2 + a3 * p3 + a4 * p4
    //   = q
    // 1 = a1 + a2 + a3 + a4
    // Da1 = [ q p2 p3 p4 ]
    // Da2 = [ p1 q p3 p4 ]
    // Da3 = [ p1 p2 q p4 ]
    // Da4 = [ p1 p2 p3 q ]
    tetrahedron_coordinates_numerator[1]=Adaptive_Signed_Volume<T_EXACT>(point,tetrahedron[2],tetrahedron[3],tetrahedron[4]);
    tetrahedron_coordinates_numerator[2]=Adaptive_Signed_Volume<T_EXACT>(tetrahedron[1],point,tetrahedron[3],tetrahedron[4]);
    tetrahedron_coordinates_numerator[3]=Adaptive_Signed_Volume<T_EXACT>(tetrahedron[1],tetrahedron[2],point,tetrahedron[4]);
    tetrahedron_coordinates_numerator[4]=Adaptive_Signed_Volume<T_EXACT>(tetrahedron[1],tetrahedron[2],tetrahedron[3],point);
    // Da  = [ p1 p2 p3 p4 ]
    // a1  = Da1 / Da
    // a2  = Da2 / Da
    // a3  = Da3 / Da
    // a4  = Da4 / Da
}

template<class T_EXACT,class T,class A1,class A2>
void Intersection_Coordinates_Helper(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
    A1& Da1,A1& Da2,A1& Da3,A2& Da)
{
    STATIC_ASSERT((IS_ADAPTIVE<A1>::value));
    STATIC_ASSERT((IS_ADAPTIVE<A2>::value));
    typedef INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE1 T_ADAPTIVE1;
    // triangle1 = [p1 p2 p3]
    // triangle2 = [q1 q2 q3]
    // triangle3 = [r1 r2 r3]
    // x = a1 * p1 + a2 * p2 + a3 * p3
    //   = b1 * q1 + b2 * q2 + b3 * q3
    //   = c1 * r1 + c2 * r2 + c3 * r3
    // 1 = a1 + a2 + a3
    //   = b1 + b2 + b3
    //   = c1 + c2 + c3
    // A11 = 1       A12 = 1       A13 = 1
    // A21 = [p1 q]  A22 = [p2 q]  A23 = [p3 q]
    // A31 = [p1 r]  A32 = [p2 r]  A33 = [p3 r]
    T_ADAPTIVE1 A21(Adaptive_Signed_Volume<T_EXACT>(triangle1[1],triangle2));
    T_ADAPTIVE1 A22(Adaptive_Signed_Volume<T_EXACT>(triangle1[2],triangle2));
    T_ADAPTIVE1 A23(Adaptive_Signed_Volume<T_EXACT>(triangle1[3],triangle2));
    T_ADAPTIVE1 A31(Adaptive_Signed_Volume<T_EXACT>(triangle1[1],triangle3));
    T_ADAPTIVE1 A32(Adaptive_Signed_Volume<T_EXACT>(triangle1[2],triangle3));
    T_ADAPTIVE1 A33(Adaptive_Signed_Volume<T_EXACT>(triangle1[3],triangle3));
    Da1=Adaptive_Determinant<void>(A22,A23,A32,A33);
    Da2=Adaptive_Determinant<void>(A23,A21,A33,A31);
    Da3=Adaptive_Determinant<void>(A21,A22,A31,A32);
    Da=Da1+Da2+Da3;
    // a1 = Da1 / Da
    // a2 = Da2 / Da
    // a3 = Da3 / Da
}

//#################################################################
// Intersection_Coordinates
//#################################################################
template<class T_EXACT,int N,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<T,N>& point,const VECTOR<VECTOR<T,N>,N+1>& simplex,VECTOR<T_ADAPTIVE,N+1>& simplex_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    Intersection_Coordinates<T_EXACT>(simplex,point,simplex_coordinates);
}

template<class T_EXACT,int N,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<VECTOR<T,N>,N+1>& simplex,const VECTOR<T,N>& point,VECTOR<T_ADAPTIVE,N+1>& simplex_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    STATIC_ASSERT((N==2 || N==3));
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,N>,N+1>,VECTOR<T,N> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE2 ADAPTIVE_TYPE2;
    ADAPTIVE_TYPE2 common_denominator(Adaptive_Signed_Volume<T_EXACT>(simplex));
    VECTOR<ADAPTIVE_TYPE1,N+1> simplex_coordinates_numerator;
    Intersection_Coordinates_Helper<T_EXACT>(simplex,point,simplex_coordinates_numerator);
    for(int i=1;i<=N+1;++i) simplex_coordinates[i]=simplex_coordinates_numerator[i]/common_denominator;
}

template<class T_EXACT,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<VECTOR<T,2>,2>& segment1,const VECTOR<VECTOR<T,2>,2>& segment2,VECTOR<T_ADAPTIVE,2>& segment1_coordinates,VECTOR<T_ADAPTIVE,2>& segment2_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    Intersection_Coordinates<T_EXACT>(segment1,segment2,segment1_coordinates);
    Intersection_Coordinates<T_EXACT>(segment2,segment1,segment2_coordinates);
}

template<class T_EXACT,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<VECTOR<T,2>,2>& segment1,const VECTOR<VECTOR<T,2>,2>& segment2,VECTOR<T_ADAPTIVE,2>& segment1_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,2>,2>,VECTOR<VECTOR<T,2>,2> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE2 ADAPTIVE_TYPE2;
    // segment1 = [p1 p2]
    // segment2 = [q1 q2]
    // x = a1 * p1 + a2 * p2
    //   = b1 * q1 + b2 * q2
    // 1 = a1 + a2
    //   = b1 + b2
    // A11 = 1       A12 = 1
    // A21 = [p1 q]  A22 = [p2 q]
    ADAPTIVE_TYPE1 Da1(Adaptive_Signed_Volume<T_EXACT>(segment1[2],segment2[1],segment2[2]));
    ADAPTIVE_TYPE1 Da2(Adaptive_Signed_Volume<T_EXACT>(segment1[1],segment2[2],segment2[1]));
    ADAPTIVE_TYPE2 Da=Da1+Da2;
    segment1_coordinates[1]=Da1/Da;
    segment1_coordinates[2]=Da2/Da;
}

template<class T_EXACT,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<VECTOR<T,3>,2>& segment,const VECTOR<VECTOR<T,3>,3>& triangle,VECTOR<T_ADAPTIVE,2>& segment_coordinates,VECTOR<T_ADAPTIVE,3>& triangle_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    Intersection_Coordinates<T_EXACT>(triangle,segment,triangle_coordinates,segment_coordinates);
}

template<class T_EXACT,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment,VECTOR<T_ADAPTIVE,3>& triangle_coordinates,VECTOR<T_ADAPTIVE,2>& segment_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,2> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE2 ADAPTIVE_TYPE2;
    // triangle = [p1 p2 p3]
    // segment  = [q1 q2]
    // x = a1 * p1 + a2 * p2 + a3 * p3
    //   = b1 * q1 + b2 * q2
    // 1 = a1 + a2 + a3
    //   = b1 + b2
    // Da1 =  [p2 p3 q1 q2]
    // Da2 = -[p1 p3 q1 q2] = [p3 p1 q1 q2]
    // Da3 =  [p1 p2 q1 q2]
    // Db1 =  [p1 p2 p3 q2]
    // Db2 = -[p1 p2 p3 q1] = [p1 p3 p2 q1]
    // Dab =  [p1 p2 p3 q2] - [p1 p2 p3 q1] = Db1 + Db2
    ADAPTIVE_TYPE1 Da1(Adaptive_Signed_Volume<T_EXACT>(triangle[2],triangle[3],segment[1],segment[2]));
    ADAPTIVE_TYPE1 Da2(Adaptive_Signed_Volume<T_EXACT>(triangle[3],triangle[1],segment[1],segment[2]));
    ADAPTIVE_TYPE1 Da3(Adaptive_Signed_Volume<T_EXACT>(triangle[1],triangle[2],segment[1],segment[2]));
    ADAPTIVE_TYPE1 Db1(Adaptive_Signed_Volume<T_EXACT>(triangle[1],triangle[2],triangle[3],segment[2]));
    ADAPTIVE_TYPE1 Db2(Adaptive_Signed_Volume<T_EXACT>(triangle[1],triangle[3],triangle[2],segment[1]));
    ADAPTIVE_TYPE2 Dab(Db1+Db2);
    triangle_coordinates[1]=Da1/Dab;
    triangle_coordinates[2]=Da2/Dab;
    triangle_coordinates[3]=Da3/Dab;
    segment_coordinates[1]=Db1/Dab;
    segment_coordinates[2]=Db2/Dab;
}

template<class T_EXACT,class T,class T_ADAPTIVE> void
Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,VECTOR<T_ADAPTIVE,3>& triangle1_coordinates,
    VECTOR<T_ADAPTIVE,3>& triangle2_coordinates,VECTOR<T_ADAPTIVE,3>& triangle3_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    Intersection_Coordinates<T_EXACT>(triangle1,triangle2,triangle3,triangle1_coordinates);
    Intersection_Coordinates<T_EXACT>(triangle2,triangle3,triangle1,triangle2_coordinates);
    Intersection_Coordinates<T_EXACT>(triangle3,triangle1,triangle2,triangle3_coordinates);
}

template<class T_EXACT,class T,class T_ADAPTIVE> void
Intersection_Coordinates(VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,VECTOR<T_ADAPTIVE,3>& triangle1_coordinates)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::ADAPTIVE_TYPE2 ADAPTIVE_TYPE2;
    ADAPTIVE_TYPE1 Da1,Da2,Da3;
    ADAPTIVE_TYPE2 Da;
    Intersection_Coordinates_Helper<T_EXACT>(triangle1,triangle2,triangle3,Da1,Da2,Da3,Da);
    triangle1_coordinates[1]=Da1/Da;
    triangle1_coordinates[2]=Da2/Da;
    triangle1_coordinates[3]=Da3/Da;
}

//#################################################################
// Intersection_Coordinates
// computes non-adaptive coordinates to an absolute tolerance
//#################################################################
template<class T_EXACT,int N,class T,class U,class V> void
Intersection_Coordinates(const VECTOR<T,N>& point,const VECTOR<VECTOR<T,N>,N+1>& simplex,VECTOR<U,N+1>& simplex_coordinates,V abs_tol)
{
    Intersection_Coordinates<T_EXACT>(simplex,point,simplex_coordinates,abs_tol);
}

template<class T_EXACT,int N,class T,class U,class V> void
Intersection_Coordinates(const VECTOR<VECTOR<T,N>,N+1>& simplex,const VECTOR<T,N>& point,VECTOR<U,N+1>& simplex_coordinates,V abs_tol)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<U>::value));
    STATIC_ASSERT((N == 2 || N == 3));
    assert(abs_tol > 0);
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,N>,N+1>,VECTOR<T,N> > ADAPTIVE_RESULT;typedef typename ADAPTIVE_RESULT::TYPE ADAPTIVE_TYPE;
    VECTOR<ADAPTIVE_TYPE,N+1> adaptive_coordinates;
    Intersection_Coordinates<T_EXACT>(simplex,point,adaptive_coordinates);
    Convert_And_Ensure_Abs_Tol(adaptive_coordinates,simplex_coordinates,abs_tol);
}

template<class T_EXACT,class T,class U,class V> void
Intersection_Coordinates(const VECTOR<VECTOR<T,2>,2>& simplex1,const VECTOR<VECTOR<T,2>,2>& simplex2,VECTOR<U,2>& simplex1_coordinates,V abs_tol)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<U>::value));
    assert(abs_tol > 0);
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,2>,2>,VECTOR<VECTOR<T,2>,2> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE ADAPTIVE_TYPE;
    VECTOR<ADAPTIVE_TYPE,2> adaptive1_coordinates;
    Intersection_Coordinates<T_EXACT>(simplex1,simplex2,adaptive1_coordinates);
    Convert_And_Ensure_Abs_Tol(adaptive1_coordinates,simplex1_coordinates,abs_tol);
}

template<class T_EXACT,int D,int N1,int N2,class T,class U,class V> void
Intersection_Coordinates(const VECTOR<VECTOR<T,D>,N1>& simplex1,const VECTOR<VECTOR<T,D>,N2>& simplex2,VECTOR<U,N1>& simplex1_coordinates,VECTOR<U,N2>& simplex2_coordinates,V abs_tol)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<U>::value));
    STATIC_ASSERT(((D == 2 && N1 == 2 && N2 == 2) || (D == 3 && ((N1 == 2 && N2 == 3) || (N1 == 3 && N2 == 2)))));
    assert(abs_tol > 0);
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,D>,N1>,VECTOR<VECTOR<T,D>,N2> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE ADAPTIVE_TYPE;
    VECTOR<ADAPTIVE_TYPE,N1> adaptive1_coordinates;
    VECTOR<ADAPTIVE_TYPE,N2> adaptive2_coordinates;
    Intersection_Coordinates<T_EXACT>(simplex1,simplex2,adaptive1_coordinates,adaptive2_coordinates);
    Convert_And_Ensure_Abs_Tol(adaptive1_coordinates,simplex1_coordinates,abs_tol);
    Convert_And_Ensure_Abs_Tol(adaptive2_coordinates,simplex2_coordinates,abs_tol);
}

template<class T_EXACT,class T,class U,class V> void
Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& simplex1,const VECTOR<VECTOR<T,3>,3>& simplex2,const VECTOR<VECTOR<T,3>,3>& simplex3,VECTOR<U,3>& simplex1_coordinates,V abs_tol)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<T>::value));
    assert(abs_tol > 0);
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE ADAPTIVE_TYPE;
    VECTOR<ADAPTIVE_TYPE,3> adaptive1_coordinates;
    Intersection_Coordinates<T_EXACT>(simplex1,simplex2,simplex3,adaptive1_coordinates);
    Convert_And_Ensure_Abs_Tol(adaptive1_coordinates,simplex1_coordinates,abs_tol);
}

template<class T_EXACT,class T,class U,class V> void
Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& simplex1,const VECTOR<VECTOR<T,3>,3>& simplex2,const VECTOR<VECTOR<T,3>,3>& simplex3,
    VECTOR<U,3>& simplex1_coordinates,VECTOR<U,3>& simplex2_coordinates,VECTOR<U,3>& simplex3_coordinates,V abs_tol)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<U>::value));
    assert(abs_tol > 0);
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE ADAPTIVE_TYPE;
    VECTOR<ADAPTIVE_TYPE,3> adaptive1_coordinates;
    VECTOR<ADAPTIVE_TYPE,3> adaptive2_coordinates;
    VECTOR<ADAPTIVE_TYPE,3> adaptive3_coordinates;
    Intersection_Coordinates<T_EXACT>(simplex1,simplex2,simplex3,adaptive1_coordinates,adaptive2_coordinates,adaptive3_coordinates);
    Convert_And_Ensure_Abs_Tol(adaptive1_coordinates,simplex1_coordinates,abs_tol);
    Convert_And_Ensure_Abs_Tol(adaptive2_coordinates,simplex2_coordinates,abs_tol);
    Convert_And_Ensure_Abs_Tol(adaptive3_coordinates,simplex3_coordinates,abs_tol);
}
//#################################################################
// Convert_And_Ensure_Abs_Tol
//#################################################################
template<class T_ADAPTIVE,class T,int N,class U>
void Convert_And_Ensure_Abs_Tol(const VECTOR<T_ADAPTIVE,N>& adaptive_coordinates,VECTOR<T,N>& simplex_coordinates,U abs_tol)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef typename T_ADAPTIVE::FP_TYPE FP_TYPE;
    for(int i=1;i<=N;++i){
        FP_TYPE estimate,error;
        adaptive_coordinates[i].Estimate_And_Error().Get(estimate,error);
        if(error>abs_tol) estimate=adaptive_coordinates[i].Refined_Estimate();
        simplex_coordinates[i]=estimate;}
}
//#####################################################################
}
using GEOMETRIC_PREDICATES_DETAIL::Intersection_Coordinates;
using GEOMETRIC_PREDICATES_DETAIL::INTERSECTION_COORDINATES_ADAPTIVE_RESULT;
}
#endif
