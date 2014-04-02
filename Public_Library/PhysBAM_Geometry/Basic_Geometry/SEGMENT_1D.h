//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_1D  
//##################################################################### 
#ifndef __SEGMENT_1D__
#define __SEGMENT_1D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class TV> class RANGE;

template<class T>
class SEGMENT_1D
{
    typedef VECTOR<T,1> TV;
public:
    VECTOR<T,1> x1,x2;

    SEGMENT_1D()
        :x1((T)0),x2((T)1)
    {}

    SEGMENT_1D(const VECTOR<T,1>& x1_input,const VECTOR<T,1>& x2_input)
        :x1(x1_input),x2(x2_input)
    {}

    template<class T_ARRAY>
    explicit SEGMENT_1D(const T_ARRAY& X_input)
        :x1(X_input(1)),x2(X_input(2))
    {
        STATIC_ASSERT(T_ARRAY::m==2);
    }

    T Length() const
    {return (x2-x1).Magnitude();}

    T Size() const
    {return Length();}

    VECTOR<T,2> Barycentric_Coordinates(const TV& location) const
    {return Barycentric_Coordinates(location,x1,x2);}

    VECTOR<T,2> Sum_Barycentric_Coordinates(const SEGMENT_1D<T>& embedded_segment) const
    {return Barycentric_Coordinates(embedded_segment.x1)+Barycentric_Coordinates(embedded_segment.x2);}

    static VECTOR<T,2> Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2)
    {T fraction=(location.x-x1.x)/(x2-x1).Magnitude();return VECTOR<T,2>((T)1-fraction,fraction);}

    template<class T_ARRAY>
    static VECTOR<T,2> Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Barycentric_Coordinates(location,X(1),X(2));}

    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights,const TV& x1,const TV& x2)
    {return weights.x*x1+weights.y*x2;}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x1,const TV& x2)
    {return weights.x*x1+(1-weights.x)*x2;}

    TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights) const
    {return Point_From_Barycentric_Coordinates(weights,x1,x2);}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,2>& weights,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Point_From_Barycentric_Coordinates(weights,X(1),X(2));}

    VECTOR<T,1> Center() const
    {return (T).5*(x1+x2);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(x1,x2);}

    bool Inside(const TV& location,const T thickness_over_2=0)
    {return x1.x<=location.x && location.x<=x2.x;}

//#####################################################################
    static T Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,2>& indices){PHYSBAM_NOT_IMPLEMENTED();}
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<SEGMENT_1D<T> >& clipped_simplices) const {PHYSBAM_NOT_IMPLEMENTED();}
    VECTOR<T,1> Closest_Point_On_Segment(const VECTOR<T,1>& point) const {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};

template<class T> std::ostream &operator<<(std::ostream &output,const SEGMENT_1D<T> &segment)
{output << segment.x1 << ", " << segment.x2;return output;}

}
#endif

