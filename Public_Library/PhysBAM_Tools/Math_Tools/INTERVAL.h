//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERVAL
//#####################################################################
#ifndef __INTERVAL__
#define __INTERVAL__

#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/ZERO.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <cfloat>
namespace PhysBAM{

template<class TV> struct IS_SCALAR_BLOCK<INTERVAL<TV> >:public IS_SCALAR_BLOCK<TV>{};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<INTERVAL<TV>,RW> {static const bool value=true;};

template<class T>
class INTERVAL
{
public:
    template<class T2> struct REBIND{typedef INTERVAL<T2> TYPE;};
    typedef T SCALAR;

    T min_corner,max_corner;

    INTERVAL()
        :min_corner(),max_corner()
    {}

    INTERVAL(const T value)
        :min_corner(value),max_corner(value)
    {}

    INTERVAL(const T min_corner,const T max_corner)
        :min_corner(min_corner),max_corner(max_corner)
    {}

    template<class T2> explicit INTERVAL(const INTERVAL<T2>& interval)
        :min_corner(T(interval.min_corner)),max_corner(T(interval.max_corner))
    {}

    static INTERVAL Unit_Box()
    {return INTERVAL((T)0,(T)1);}

    static INTERVAL Zero_Box()
    {return INTERVAL();}

    static INTERVAL Empty_Box()
    {return INTERVAL((T)FLT_MAX,-(T)FLT_MAX);}

    static INTERVAL Full_Box()
    {return INTERVAL(-(T)FLT_MAX,(T)FLT_MAX);}

    bool Empty() const
    {return min_corner>max_corner;}

    bool operator==(const INTERVAL& r) const
    {return min_corner==r.min_corner && max_corner==r.max_corner;}

    bool operator!=(const INTERVAL& r) const
    {return !(*this==r);}

    INTERVAL operator-() const
    {return INTERVAL(-max_corner,-min_corner);}

    INTERVAL& operator+=(const INTERVAL& r)
    {min_corner+=r.min_corner;max_corner+=r.max_corner;return *this;}

    INTERVAL& operator-=(const INTERVAL& r)
    {min_corner-=r.max_corner;max_corner-=r.min_corner;return *this;}

    INTERVAL operator+(const INTERVAL& r) const
    {return INTERVAL(min_corner+r.min_corner,max_corner+r.max_corner);}

    INTERVAL operator-(const INTERVAL& r) const
    {return INTERVAL(min_corner-r.max_corner,max_corner-r.min_corner);}

    INTERVAL operator*(const T a) const
    {return a>=0?INTERVAL(min_corner*a,max_corner*a):INTERVAL(max_corner*a,min_corner*a);}

    INTERVAL& operator*=(const T a)
    {return *this=*this*a;}

    INTERVAL operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    INTERVAL& operator/=(const T a)
    {return *this=*this/a;}

    T Center() const
    {return (T).5*(min_corner+max_corner);}

    T Size() const
    {return max_corner-min_corner;}

    void Enlarge_To_Include_Point(const T& point)
    {min_corner=min(min_corner,point);max_corner=max(max_corner,point);}

    void Enlarge_Nonempty_Box_To_Include_Point(const T& point)
    {assert(!Empty());if(point<min_corner) min_corner=point;else if(point>max_corner) max_corner=point;}

    void Enlarge_Nonempty_Box_To_Include_Points(const T& p1,const T& p2)
    {Enlarge_Nonempty_Box_To_Include_Point(p1);Enlarge_Nonempty_Box_To_Include_Point(p2);}

    void Enlarge_Nonempty_Box_To_Include_Points(const T& p1,const T& p2,const T& p3)
    {Enlarge_Nonempty_Box_To_Include_Point(p1);Enlarge_Nonempty_Box_To_Include_Point(p2);Enlarge_Nonempty_Box_To_Include_Point(p3);}

    template<class T_ARRAY>
    void Enlarge_Nonempty_Box_To_Include_Points(const T_ARRAY& points)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=1;i<=points.Size();i++) Enlarge_Nonempty_Box_To_Include_Point(points(i));}

    void Enlarge_To_Include_Box(const INTERVAL& interval)
    {min_corner=min(min_corner,interval.min_corner);max_corner=max(max_corner,interval.max_corner);}

    void Change_Size(const T delta)
    {min_corner-=delta;max_corner+=delta;}

    INTERVAL Thickened(const T thickness_over_two) const
    {return INTERVAL(min_corner-thickness_over_two,max_corner+thickness_over_two);}

    static INTERVAL Combine(const INTERVAL& box1,const INTERVAL& box2)
    {return INTERVAL(min(box1.min_corner,box2.min_corner),max(box1.max_corner,box2.max_corner));}

    static INTERVAL Intersect(const INTERVAL& box1,const INTERVAL& box2)
    {return INTERVAL(max(box1.min_corner,box2.min_corner),min(box1.max_corner,box2.max_corner));}

    void Scale_About_Center(const T factor)
    {T center=(T).5*(min_corner+max_corner),length_over_two=factor*(T).5*(max_corner-min_corner);min_corner=center-length_over_two;max_corner=center+length_over_two;}

    bool Lazy_Inside(const T& location) const
    {return min_corner<=location && location<=max_corner;}

    bool Lazy_Inside_Half_Open(const T& location) const
    {return min_corner<=location && location<max_corner;}

    bool Inside(const T& location,const T thickness_over_two) const
    {return Thickened(-thickness_over_two).Lazy_Inside(location);}

    bool Inside(const T& location,const ZERO thickness_over_two) const
    {return Lazy_Inside(location);}

    bool Lazy_Outside(const T& location) const
    {return !Lazy_Inside(location);}

    bool Outside(const T& location,const T thickness_over_two) const
    {return Thickened(thickness_over_two).Lazy_Outside(location);}

    bool Outside(const T& location,const ZERO thickness_over_two) const
    {return Lazy_Outside(location);}

    bool Boundary(const T& location,const T thickness_over_two) const
    {bool strict_inside=min_corner+thickness_over_two<location && location<max_corner-thickness_over_two;
    return !strict_inside && !Outside(location,thickness_over_two);}

    T Clamp(const T& location) const
    {return clamp(location,min_corner,max_corner);}

    bool Contains(const INTERVAL& interval) const
    {return min_corner<=interval.min_corner && interval.max_corner<=max_corner;}

    bool Lazy_Intersection(const INTERVAL& interval) const
    {return min_corner<=interval.max_corner && interval.min_corner<=max_corner;}

    bool Intersection(const INTERVAL& interval,const T thickness_over_two) const
    {return Thickened(thickness_over_two).Lazy_Intersection(interval);}

    bool Intersection(const INTERVAL& interval,const ZERO thickness_over_two) const
    {return Lazy_Intersection(interval);}

    bool Intersection(const INTERVAL& interval) const
    {return Lazy_Intersection(interval);}

    static INTERVAL Bounding_Box(const T& p1,const T& p2)
    {INTERVAL interval(p1);interval.Enlarge_Nonempty_Box_To_Include_Point(p2);return interval;}

    static INTERVAL Bounding_Box(const T& p1,const T& p2,const T& p3)
    {INTERVAL interval(p1);interval.Enlarge_Nonempty_Box_To_Include_Points(p2,p3);return interval;}

    static INTERVAL Bounding_Box(const T& p1,const T& p2,const T& p3,const T& p4)
    {INTERVAL interval(p1);interval.Enlarge_Nonempty_Box_To_Include_Points(p2,p3,p4);return interval;}

    template<class T_ARRAY>
    static INTERVAL Bounding_Box(const T_ARRAY& points)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    if(!points.Size()) return Empty_Box();
    INTERVAL interval(points(1));for(int i=2;i<=points.Size();i++) interval.Enlarge_Nonempty_Box_To_Include_Point(points(i));return interval;}

    T Signed_Distance(const T& X) const
    {return abs(X-Center())-Size();}

//#####################################################################
};

template<class T>
inline INTERVAL<T> operator+(const T& a,const INTERVAL<T>& b)
{return INTERVAL<T>(a+b.min_corner,a+b.max_corner);}

template<class T>
inline INTERVAL<T> operator-(const T& a,const INTERVAL<T>& b)
{return INTERVAL<T>(a-b.max_corner,a-b.min_corner);}

template<class T> inline INTERVAL<T> operator*(const typename T::SCALAR a,const INTERVAL<T>& interval)
{return interval*a;}
}
#endif
