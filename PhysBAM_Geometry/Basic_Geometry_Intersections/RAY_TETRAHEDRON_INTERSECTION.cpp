//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TETRAHEDRON_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const TETRAHEDRON<T>& tetrahedron, const T thickness)
{
    bool intersection=false;
    if(INTERSECTION::Intersects(ray,tetrahedron.triangle1,thickness)){intersection=true;ray.aggregate_id=1;}
    if(INTERSECTION::Intersects(ray,tetrahedron.triangle2,thickness)){intersection=true;ray.aggregate_id=2;}
    if(INTERSECTION::Intersects(ray,tetrahedron.triangle3,thickness)){intersection=true;ray.aggregate_id=3;}
    if(INTERSECTION::Intersects(ray,tetrahedron.triangle4,thickness)){intersection=true;ray.aggregate_id=4;}
    if(intersection) return true;
    else return false;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const TETRAHEDRON<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const TETRAHEDRON<double>&,const double);
#endif
};
};
