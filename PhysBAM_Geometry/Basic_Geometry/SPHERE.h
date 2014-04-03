//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE
//##################################################################### 
#ifndef __SPHERE__
#define __SPHERE__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class TV>
class SPHERE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:
    typedef TV VECTOR_T;

    TV center;
    T radius;

    SPHERE()
        :radius(1)
    {}

    SPHERE(const TV& center,const T radius)
        :center(center),radius(radius)
    {}

    TV Normal(const TV& location) const
    {return (location-center).Normalized();}

    bool Inside(const TV& location,const T thickness_over_two) const
    {return (location-center).Magnitude_Squared() <= sqr(radius-thickness_over_two);}

    bool Lazy_Inside(const TV& location) const
    {return (location-center).Magnitude_Squared() <= sqr(radius);}
    
    bool Outside(const TV& location,const T thickness_over_two) const
    {return (location-center).Magnitude_Squared() >= sqr(radius+thickness_over_two);}

    bool Lazy_Outside(const TV& location) const
    {return (location-center).Magnitude_Squared() >= sqr(radius);}

    bool Boundary(const TV& location,const T thickness_over_two) const
    {return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);}

    TV Surface(const TV& location) const  
    {return radius*(location-center).Normalized()+center;}

    T Signed_Distance(const TV& location) const
    {return (location-center).Magnitude()-radius;}

    T Circular_Segment_Area(const T h) const
    {STATIC_ASSERT(d==2);return sqr(radius)*acos((radius-h)/radius)-(radius-h)*sqrt(2*radius*h-sqr(h));}

    T Size() const
    {return (T)unit_sphere_size<d>::value*pow<d>(radius);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(center).Thickened(radius);}

    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,d-1>::All_Ones_Vector()/radius;}

//#####################################################################
    void Sector_Volumes(const TV& origin,T volumes[1<<d],const T thickness_over_two=0) const;
    T Octant_Volume(const VECTOR<T,3>& min_corner) const;
    static std::string Name();
//#####################################################################
};   
}
#endif
