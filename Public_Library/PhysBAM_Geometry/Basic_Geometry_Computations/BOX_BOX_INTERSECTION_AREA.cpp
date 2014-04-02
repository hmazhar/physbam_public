//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersection_Area
//#####################################################################
template<class TV> typename TV::SCALAR Intersection_Area(const RANGE<TV>& box1, const RANGE<TV>& box2)
{
    return box1.Intersection_Area(box2);
}
//#####################################################################
template float Intersection_Area(const RANGE<VECTOR<float,1> >&,const RANGE<VECTOR<float,1> >&);
template float Intersection_Area(const RANGE<VECTOR<float,2> >&,const RANGE<VECTOR<float,2> >&);
template float Intersection_Area(const RANGE<VECTOR<float,3> >&,const RANGE<VECTOR<float,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template double Intersection_Area(const RANGE<VECTOR<double,1> >&,const RANGE<VECTOR<double,1> >&);
template double Intersection_Area(const RANGE<VECTOR<double,2> >&,const RANGE<VECTOR<double,2> >&);
template double Intersection_Area(const RANGE<VECTOR<double,3> >&,const RANGE<VECTOR<double,3> >&);
#endif
};
};
