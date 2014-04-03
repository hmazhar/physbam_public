//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CYLINDER
//##################################################################### 
#ifndef __CYLINDER__
#define __CYLINDER__

#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{

template<class T>
class CYLINDER
{
    typedef VECTOR<T,3> TV;
public:
    typedef TV VECTOR_T;

    T height,radius;
    PLANE<T> plane1,plane2; // plane2 is height units behind circle

    CYLINDER()
        :height(1),radius(1),plane2(plane1.x1-height*plane1.normal,-plane1.normal)
    {}

    CYLINDER(const TV& point1,const TV& point2,const T radius)
        :radius(radius)
    {
        Set_Endpoints(point1,point2);
    }

    void Set_Height(const T height_input)
    {height=height_input;plane2.x1=plane1.x1-height*plane1.normal;}

    void Set_Endpoints(const TV& point1,const TV& point2)
    {TV height_vector=point1-point2;height=height_vector.Normalize();
    plane1.x1=point1;plane1.normal=height_vector;
    plane2.x1=point2;plane2.normal=-plane1.normal;}

    VECTOR<T,2> Principal_Curvatures(const TV& X) const
    {TV v=X-plane1.x1;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane_distance=max(plane1_distance,-height-plane1_distance);
    T cylinder_distance=(v-plane1_distance*plane1.normal).Magnitude()-radius;
    if(abs(plane_distance)<abs(cylinder_distance)) return VECTOR<T,2>();
    return VECTOR<T,2>(-1/radius,0);}

    static std::string Name() 
    {return "CYLINDER<T>";}

//#####################################################################
    TV Normal(const TV& location) const;
    TV Normal(const TV& location,const int aggregate) const;
    bool Inside(const TV& location,const T thickness_over_two) const;
    bool Lazy_Inside(const TV& location) const;
    bool Outside(const TV& location,const T thickness_over_two) const;
    bool Lazy_Outside(const TV& location) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    TV Surface(const TV& location) const;
    T Signed_Distance(const TV& location) const;
    RANGE<TV> Bounding_Box() const;
//#####################################################################
};   
}
#endif
