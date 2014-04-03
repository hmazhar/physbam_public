//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERSECTS_AND_INTERSECTION_COORDINATES__
#define __INTERSECTS_AND_INTERSECTION_COORDINATES__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OBJECT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_SIGNED_VOLUME.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersection_Coordinates.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects_And_Intersection_Coordinates_Forward.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects_Forward.h>

namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
//#################################################################
// Intersects_And_Intersection_Coordinates_Helper
//#################################################################
template<class T_EXACT,class T,class T_ADAPTIVE>
bool Intersects_And_Intersection_Coordinates_Helper(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
    VECTOR<T_ADAPTIVE,3>& triangle1_coordinates,bool* is_degenerate_p)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::TYPE2 ADAPTIVE_TYPE2;
    if(is_degenerate_p!=0) *is_degenerate_p=false;
    ADAPTIVE_TYPE1 Dx1,Dx2,Dx3;
    ADAPTIVE_TYPE2 Dx;
    Intersection_Coordinates_Helper<T_EXACT>(triangle1,triangle2,triangle3,Dx1,Dx2,Dx3,Dx);
    if(Dx1.Sign()!=0 && Dx2.Sign()!=0 && Dx3.Sign()!=0 && Dx.Sign()!=0){
        // barycentric coordinates for intersection are non-degenerate, so see if they are all the same sign
        if(Dx1.Sign()!=Dx2.Sign() || Dx1.Sign()!=Dx3.Sign()) return false;
        triangle1_coordinates[1]=Dx1/Dx;triangle1_coordinates[2]=Dx2/Dx;triangle1_coordinates[3]=Dx3/Dx;
        return true;}
    // barycentric coordinates are degenerate, but perhaps we can determine that the 3 triangles couldn't intersect anyway...
    bool is_degenerate;
    if(!Intersects<T_EXACT>(triangle1,triangle2,&is_degenerate) && !is_degenerate) return false;
    if(!Intersects<T_EXACT>(triangle1,triangle3,&is_degenerate) && !is_degenerate) return false;
    if(!Intersects<T_EXACT>(triangle2,triangle3,&is_degenerate) && !is_degenerate) return false;
    if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();*is_degenerate_p=true;
    return true;
}
//#################################################################
// Intersects_And_Intersection_Coordinates
//#################################################################
// If there is an intersection, computes adaptive barycentric coordinates of the intersection.
template<class T_EXACT,class T,class T_ADAPTIVE>
bool Intersects_And_Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
    VECTOR<T_ADAPTIVE,3>& triangle1_coordinates,VECTOR<T_ADAPTIVE,3>& triangle2_coordinates,VECTOR<T_ADAPTIVE,3>& triangle3_coordinates,bool* is_degenerate_p)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    typedef INTERSECTION_COORDINATES_HELPER_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE1 ADAPTIVE_TYPE1;
    typedef typename ADAPTIVE_RESULT::TYPE2 ADAPTIVE_TYPE2;
    ADAPTIVE_TYPE1 Dx1,Dx2,Dx3;
    ADAPTIVE_TYPE2 Dx;
    if(!Intersects_And_Intersection_Coordinates_Helper<T_EXACT>(triangle1,triangle2,triangle3,triangle1_coordinates,is_degenerate_p)) return false;
    if(is_degenerate_p!=0&&*is_degenerate_p) return true; // degenerate, so no need to continue
    if(!Intersects_And_Intersection_Coordinates_Helper<T_EXACT>(triangle2,triangle3,triangle1,triangle2_coordinates,is_degenerate_p)) return false;
    if(is_degenerate_p!=0&&*is_degenerate_p) return true; // degenerate, so no need to continue
    return Intersects_And_Intersection_Coordinates_Helper<T_EXACT>(triangle3,triangle1,triangle2,triangle3_coordinates,is_degenerate_p);
}
//#################################################################
// Intersects_And_Intersection_Coordinates
//#################################################################
// If there is an intersection, computes non-adaptive barycentric coordinates of the intersection to an absolute tolerance.
template<class T_EXACT,class T,class U,class V>
bool Intersects_And_Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
    VECTOR<U,3>& triangle1_coordinates,VECTOR<U,3>& triangle2_coordinates,VECTOR<U,3>& triangle3_coordinates,V abs_tol,bool* is_degenerate_p)
{
    STATIC_ASSERT((IS_NOT_ADAPTIVE<U>::value));
    typedef INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> > ADAPTIVE_RESULT;
    typedef typename ADAPTIVE_RESULT::TYPE ADAPTIVE_TYPE;
    assert(abs_tol>0);
    VECTOR<ADAPTIVE_TYPE,3> adaptive1_coordinates,adaptive2_coordinates,adaptive3_coordinates;
    bool intersects=Intersects_And_Intersection_Coordinates<T_EXACT>(triangle1,triangle2,triangle3,adaptive1_coordinates,adaptive2_coordinates,adaptive3_coordinates,is_degenerate_p);
    if(is_degenerate_p!=0&&*is_degenerate_p) return true; // degenerate, so no need to continue
    if(intersects){
        Convert_And_Ensure_Abs_Tol(adaptive1_coordinates,triangle1_coordinates,abs_tol);
        Convert_And_Ensure_Abs_Tol(adaptive2_coordinates,triangle2_coordinates,abs_tol);
        Convert_And_Ensure_Abs_Tol(adaptive3_coordinates,triangle3_coordinates,abs_tol);}
    return intersects;
}
//#################################################################
}
using GEOMETRIC_PREDICATES_DETAIL::Intersects_And_Intersection_Coordinates;
}
#endif
