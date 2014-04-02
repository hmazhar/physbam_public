//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TORUS
//##################################################################### 
#ifndef __TORUS__
#define __TORUS__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <string>
namespace PhysBAM{

template<class T>
class TORUS
{
    typedef VECTOR<T,3> TV;
public:
    typedef TV VECTOR_T;
    TV center,axis;
    T inner_radius,outer_radius;

    TORUS()
        :axis(TV(0,0,1)),inner_radius(1),outer_radius(2)
    {}

    TORUS(const TV& center_input,const TV& axis_input,const T inner_radius_input,const T outer_radius_input)
        :center(center_input),axis(axis_input.Normalized()),inner_radius(inner_radius_input),outer_radius(outer_radius_input)
    {}

    T Signed_Distance(const TV& X) const
    {TV x=X-center;T axial=Dot_Product(x,axis),radial=(x-axial*axis).Magnitude();return sqrt(sqr(radial-outer_radius)+sqr(axial))-inner_radius;}

    RANGE<TV> Bounding_Box() const
    {T x=sqr(axis.x),y=sqr(axis.y),z=sqr(axis.z);TV corner_offset=sqrt(TV(y+z,x+z,x+y))*outer_radius+inner_radius;
    return RANGE<TV>(center-corner_offset,center+corner_offset);}

    TV Normal(const TV& X) const
    {TV x=X-center;T axial=Dot_Product(x,axis);TV inplane=x-axial*axis;VECTOR<T,2> plane_normal(axial,inplane.Normalize()-outer_radius);plane_normal.Normalize();
    return plane_normal.x*axis+plane_normal.y*inplane;}

    VECTOR<T,2> Principal_Curvatures(const TV& X) const
    {T inplane=(X-TV::Dot_Product(X,axis)*axis).Magnitude();return VECTOR<T,2>(-1/inner_radius,-(inplane-outer_radius)/(inplane*inner_radius));}

    static std::string Name()
    {return "TORUS<T>";}

//#####################################################################
};   
}
#endif
