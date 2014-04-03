//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERSECTS_AND_INTERSECTION_COORDINATES_FORWARD__
#define __INTERSECTS_AND_INTERSECTION_COORDINATES_FORWARD__
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
//#################################################################
// Intersects_And_Intersection_Coordinates
//#################################################################
template<class EXACT_TYPE,class T,class A>
bool Intersects_And_Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
    VECTOR<A,3>& triangle1_coordinates,VECTOR<A,3>& triangle2_coordinates,VECTOR<A,3>& triangle3_coordinates,bool* is_degenerate_p=0);
template<class EXACT_TYPE,class T,class U,class V>
bool Intersects_And_Intersection_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
    VECTOR<U,3>& triangle1_coordinates,VECTOR<U,3>& triangle2_coordinates,VECTOR<U,3>& triangle3_coordinates,V abs_tol,bool* is_degenerate_p=0);
//#################################################################
}
}
#endif
