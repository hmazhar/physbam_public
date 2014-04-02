//#####################################################################
// Copyright 2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDED_HORIZONTAL_PLANE
//#####################################################################
#ifndef __BOUNDED_HORIZONTAL_PLANE__
#define __BOUNDED_HORIZONTAL_PLANE__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class TV>
class BOUNDED_HORIZONTAL_PLANE
{
    typedef typename TV::SCALAR T;
public:
    typedef TV VECTOR_T; 

    T half_width;

    BOUNDED_HORIZONTAL_PLANE(const T half_width=1)
        :half_width(half_width)
    {}

    TV Normal(const TV& location) const
    {TV normal;normal.y=1;return normal;}

    T Signed_Distance(const TV& location) const
    {return location.y;}

    bool Inside(const TV& location,const T thickness_over_two) const
    {return location.y<=-thickness_over_two;}

    bool Lazy_Inside(const TV& location) const
    {return location.y<=0;}

    bool Outside(const TV& location,const T thickness_over_two) const
    {return location.y>=thickness_over_two;}

    bool Lazy_Outside(const TV& location) const
    {return location.y>=0;}

    bool Boundary(const TV& location,const T thickness_over_two) const
    {return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);}

    TV Surface(const TV& location) const
    {TV projected=location;projected.y=0;return projected;}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(-half_width*TV::All_Ones_Vector(),half_width*TV::All_Ones_Vector());}

    VECTOR<T,TV::dimension-1> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,TV::dimension-1>();}

    static std::string Name()
    {return STRING_UTILITIES::string_sprintf("BOUNDED_HORIZONTAL_PLANE<VECTOR<T,%d> >",TV::m);}
//#####################################################################
};
}
#endif
