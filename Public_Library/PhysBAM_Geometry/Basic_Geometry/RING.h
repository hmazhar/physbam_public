//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RING
//##################################################################### 
#ifndef __RING__
#define __RING__

#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{

template<class T>
class RING
{
    typedef VECTOR<T,3> TV;
public:
    typedef TV VECTOR_T;

    T height,outer_radius,inner_radius;
    PLANE<T> plane1,plane2; // plane2 is height units behind circle

    RING()
        :height(1),outer_radius(2),inner_radius(1),plane2(plane1.x1-height*plane1.normal,-plane1.normal)
    {}

    RING(const TV& X1,const TV& X2,const T outer_radius,const T inner_radius)
        :outer_radius(outer_radius),inner_radius(inner_radius)
    {
        Set_Endpoints(X1,X2);
    }

    void Set_Height(const T height_input)
    {height=height_input;plane2.x1=plane1.x1-height*plane1.normal;}

    void Set_Endpoints(const TV& X1,const TV& X2)
    {TV plane_normal=X1-X2;height=plane_normal.Normalize();
    plane1.x1=X1;plane1.normal=plane_normal;
    plane2.x1=X2;plane2.normal=-plane_normal;}

//#####################################################################
    RANGE<TV> Bounding_Box() const;
    T Signed_Distance(const TV& X) const;
    TV Surface(const TV& X) const;
    TV Normal(const TV& X) const;
    TV Normal(const TV& X,const int aggregate) const;
    VECTOR<T,2> Principal_Curvatures(const TV& X) const;
    bool Lazy_Inside(const TV& X) const;
    bool Lazy_Outside(const TV& X) const;
    bool Inside(const TV& X,const T thickness_over_two=0) const;
    bool Outside(const TV& X,const T thickness_over_two=0) const;
    bool Boundary(const TV& X,const T thickness_over_two) const;
    static std::string Name();
//#####################################################################
};   
}
#endif
