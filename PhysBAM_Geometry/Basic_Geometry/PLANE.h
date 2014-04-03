//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANE
//#####################################################################
#ifndef __PLANE__
#define __PLANE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T>
class PLANE
{
    typedef VECTOR<T,3> TV;
public:
    typedef TV VECTOR_T;

    TV normal;
    TV x1; // point on the plane

    PLANE()
        :normal(0,1,0),x1(0,0,0)
    {}

    PLANE(const TV& normal_input,const TV& x1_input)
        :normal(normal_input),x1(x1_input)
    {}

    PLANE(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {
        Specify_Three_Points(x1_input,x2_input,x3_input);
    }

    void Specify_Three_Points(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {normal=Normal(x1_input,x2_input,x3_input);x1=x1_input;}

    static TV Normal(const TV& x1,const TV& x2,const TV& x3)
    {return TV::Cross_Product(x2-x1,x3-x1).Normalized();}

    template<class T_ARRAY>
    static TV Normal(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Normal(X(1),X(2),X(3));}

    TV Normal() const
    {return normal;}

    TV Normal(const TV& X) const
    {return normal;}

    static TV Normal_Direction(const TV& x1,const TV& x2,const TV& x3)
    {return TV::Cross_Product(x2-x1,x3-x1);} // can have any magnitude

    template<class T_ARRAY>
    static TV Normal_Direction(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Normal_Direction(X(1),X(2),X(3));}

    T Signed_Distance(const TV& location) const
    {return TV::Dot_Product(normal,location-x1);}

    // inside is the half space behind the normal
    bool Inside(const TV& location,const T thickness_over_two) const
    {return Signed_Distance(location)<=-thickness_over_two;}

    bool Lazy_Inside(const TV& location) const
    {return Signed_Distance(location)<=0;}

    bool Outside(const TV& location,const T thickness_over_two) const
    {return !Inside(location,-thickness_over_two);}

    bool Lazy_Outside(const TV& location) const
    {return !Lazy_Inside(location);}

    bool Boundary(const TV& location,const T thickness_over_two) const
    {return abs(Signed_Distance(location))<thickness_over_two;}

    // closest point on the surface
    TV Surface(const TV& location) const
    {return location-Signed_Distance(location)*normal;}

    bool Segment_Intersection(const TV& endpoint1,const TV& endpoint2,T& interpolation_fraction) const
    {return Segment_Plane_Intersection(endpoint1,endpoint2,interpolation_fraction);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Full_Box();}

    VECTOR<T,2> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,2>();}

    static std::string Name()
    {return "PLANE<T>";}

//#####################################################################
    bool Segment_Plane_Intersection(const TV& endpoint1,const TV& endpoint2,T& interpolation_fraction) const;
//#####################################################################
};
}
#endif
