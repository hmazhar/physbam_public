//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Neil Molino, Igor Neverov, Duc Nguyen, Avi Robinson-Mosher,
//     Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE
//#####################################################################
#ifndef __RANGE__
#define __RANGE__

#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/ZERO.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <cassert>
#include <cfloat>
#include <limits>
namespace PhysBAM{

template<class TV> class RANGE;
template<class TV> struct IS_SCALAR_BLOCK<RANGE<TV> >:public IS_SCALAR_BLOCK<TV>{};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<RANGE<TV>,RW> {static const bool value=false;}; // required since memory format differs from disk format

template<class TV>
class RANGE
{
    typedef typename TV::SCALAR T;
public:
    template<class T2> struct REBIND{typedef RANGE<T2> TYPE;};
    typedef T SCALAR;
    typedef TV VECTOR_T;
    enum WORKAROUND {d=TV::dimension};

    TV min_corner,max_corner;

    // Default to an empty box.
    RANGE()
        :min_corner(TV::Constant_Vector(std::numeric_limits<T>::max())),max_corner(-TV::Constant_Vector(std::numeric_limits<T>::max()))
    {}

    RANGE(const T xmin,const T xmax)
        :min_corner(xmin),max_corner(xmax)
    {
        STATIC_ASSERT(d==1);
    }

    RANGE(const T xmin,const T xmax,const T ymin,const T ymax)
        :min_corner(xmin,ymin),max_corner(xmax,ymax)
    {
        STATIC_ASSERT(d==2);
    }

    RANGE(const T xmin,const T xmax,const T ymin,const T ymax,const T zmin,const T zmax)
        :min_corner(xmin,ymin,zmin),max_corner(xmax,ymax,zmax)
    {
        STATIC_ASSERT(d==3);
    }

    RANGE(const TV& minimum_corner,const TV& maximum_corner)
        :min_corner(minimum_corner),max_corner(maximum_corner)
    {}

    template<class T2> explicit RANGE(const RANGE<T2>& box)
        :min_corner(TV(box.min_corner)),max_corner(TV(box.max_corner))
    {}

    RANGE(const TV& point)
        :min_corner(point),max_corner(point)
    {}

    static RANGE<TV> Unit_Box()
    {return RANGE<TV>(TV(),TV::All_Ones_Vector());}

    static RANGE<TV> Centered_Box()
    {return RANGE<TV>(-TV::All_Ones_Vector(),TV::All_Ones_Vector());}

    static RANGE<TV> Zero_Box()
    {return RANGE<TV>(TV(),TV());}

    static RANGE<TV> Empty_Box()
    {return RANGE<TV>();}

    static RANGE<TV> Full_Box()
    {return RANGE<TV>(TV::Constant_Vector(-std::numeric_limits<T>::max()),TV::Constant_Vector(std::numeric_limits<T>::max()));}

    bool Empty() const
    {return !min_corner.All_Less_Equal(max_corner);}

    bool operator==(const RANGE<TV>& r) const
    {return min_corner==r.min_corner && max_corner==r.max_corner;}

    bool operator!=(const RANGE<TV>& r) const
    {return !(*this==r);}

    RANGE<TV> operator-() const
    {return RANGE<TV>(-max_corner,-min_corner);}

    RANGE<TV>& operator+=(const RANGE<TV>& r)
    {min_corner+=r.min_corner;max_corner+=r.max_corner;return *this;}

    RANGE<TV>& operator-=(const RANGE<TV>& r)
    {min_corner-=r.max_corner;max_corner-=r.min_corner;return *this;}

    RANGE<TV> operator+(const RANGE<TV>& r) const
    {return RANGE<TV>(min_corner+r.min_corner,max_corner+r.max_corner);}

    RANGE<TV> operator-(const RANGE<TV>& r) const
    {return RANGE<TV>(min_corner-r.max_corner,max_corner-r.min_corner);}

    RANGE<TV> operator*(const T a) const
    {return a>=0?RANGE<TV>(min_corner*a,max_corner*a):RANGE<TV>(max_corner*a,min_corner*a);}

    RANGE<TV>& operator*=(const T a)
    {return *this=*this*a;}

    RANGE<TV> operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    RANGE<TV>& operator/=(const T a)
    {return *this=*this/a;}

    TV Edge_Lengths() const
    {return max_corner-min_corner;}

    TV Center() const
    {return (min_corner+max_corner)/2;} // Be nice to T=int

    TV Minimum_Corner() const
    {return min_corner;}

    TV Maximum_Corner() const
    {return max_corner;}

    void Corners(ARRAY<TV,VECTOR<int,1> >& corners) const
    {STATIC_ASSERT(d==1);corners.Resize(0,1);corners(0)=min_corner;corners(1)=max_corner;}

    void Corners(ARRAY<TV,VECTOR<int,2> >& corners) const
    {STATIC_ASSERT(d==2);corners.Resize(0,1,0,1);
    for(int i=0;i<=1;i++) for(int j=0;j<=1;j++)
        corners(i,j)=TV(i?max_corner.x:min_corner.x,j?max_corner.y:min_corner.y);}

    void Corners(ARRAY<TV,VECTOR<int,3> >& corners) const
    {STATIC_ASSERT(d==3);corners.Resize(0,1,0,1,0,1);
    for(int i=0;i<=1;i++) for(int j=0;j<=1;j++) for(int k=0;k<=1;k++)
        corners(i,j,k)=TV(i?max_corner.x:min_corner.x,j?max_corner.y:min_corner.y,k?max_corner.z:min_corner.z);}

    T Size() const // assumes nonnegative entries
    {return Edge_Lengths().Product();}

    T Robust_Size() const
    {return Empty()?(T)0:Size();}

    T Surface_Area() const
    {STATIC_ASSERT(d==3);VECTOR<T,3> size(Edge_Lengths());return 2*(size.x*(size.y+size.z)+size.y*size.z);}

    void Reset_Bounds(const TV& point)
    {min_corner=max_corner=point;}

    void Enlarge_To_Include_Point(const TV& point)
    {min_corner=TV::Componentwise_Min(min_corner,point);max_corner=TV::Componentwise_Max(max_corner,point);}

    void Enlarge_Nonempty_Box_To_Include_Point(const TV& point)
    {assert(!Empty());for(int i=1;i<=d;i++) if(point(i)<min_corner(i)) min_corner(i)=point(i);else if(point(i)>max_corner(i)) max_corner(i)=point(i);}

    void Enlarge_Nonempty_Box_To_Include_Points(const TV& p1,const TV& p2)
    {Enlarge_Nonempty_Box_To_Include_Point(p1);Enlarge_Nonempty_Box_To_Include_Point(p2);}

    void Enlarge_Nonempty_Box_To_Include_Points(const TV& p1,const TV& p2,const TV& p3)
    {Enlarge_Nonempty_Box_To_Include_Point(p1);Enlarge_Nonempty_Box_To_Include_Point(p2);Enlarge_Nonempty_Box_To_Include_Point(p3);}

    template<class T_ARRAY>
    void Enlarge_Nonempty_Box_To_Include_Points(const T_ARRAY& points)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    for(int i=1;i<=points.Size();i++) Enlarge_Nonempty_Box_To_Include_Point(points(i));}

    void Enlarge_To_Include_Box(const RANGE<TV>& box)
    {min_corner=TV::Componentwise_Min(min_corner,box.min_corner);max_corner=TV::Componentwise_Max(max_corner,box.max_corner);}

    void Change_Size(const T delta)
    {min_corner-=delta;max_corner+=delta;}

    void Change_Size(const TV& delta)
    {min_corner-=delta;max_corner+=delta;}

    RANGE<TV> Thickened(const T thickness_over_two) const
    {return RANGE<TV>(min_corner-thickness_over_two,max_corner+thickness_over_two);}

    static RANGE<TV> Combine(const RANGE<TV>& box1,const RANGE<TV>& box2)
    {return RANGE<TV>(TV::Componentwise_Min(box1.min_corner,box2.min_corner),TV::Componentwise_Max(box1.max_corner,box2.max_corner));}

    static RANGE<TV> Intersect(const RANGE<TV>& box1,const RANGE<TV>& box2) // assumes nonnegative entries
    {return RANGE<TV>(TV::Componentwise_Max(box1.min_corner,box2.min_corner),TV::Componentwise_Min(box1.max_corner,box2.max_corner));}

    void Scale_About_Center(const T factor)
    {TV center=(T).5*(min_corner+max_corner),length_over_two=factor*(T).5*(max_corner-min_corner);min_corner=center-length_over_two;max_corner=center+length_over_two;}

    void Scale_About_Center(const TV factor)
    {TV center=(T).5*(min_corner+max_corner),length_over_two=factor*(T).5*(max_corner-min_corner);min_corner=center-length_over_two;max_corner=center+length_over_two;}

    void Scale_About_Center(const T x_factor,const T y_factor)
    {STATIC_ASSERT(d==2);Scale_About_Center(TV(x_factor,y_factor));}

    void Scale_About_Center(const T x_factor,const T y_factor,const T z_factor)
    {STATIC_ASSERT(d==3);Scale_About_Center(TV(x_factor,y_factor,z_factor));}

    bool Lazy_Inside(const TV& location) const
    {return location.All_Greater_Equal(min_corner) && location.All_Less_Equal(max_corner);}

    bool Lazy_Inside_Half_Open(const TV& location) const
    {return location.All_Greater_Equal(min_corner) && location.All_Less(max_corner);}

    bool Inside(const TV& location,const T thickness_over_two) const
    {return Thickened(-thickness_over_two).Lazy_Inside(location);}

    bool Inside(const TV& location,const ZERO thickness_over_two) const
    {return Lazy_Inside(location);}

    bool Lazy_Outside(const TV& location) const
    {return !Lazy_Inside(location);}

    bool Outside(const TV& location,const T thickness_over_two) const
    {return Thickened(thickness_over_two).Lazy_Outside(location);}

    bool Outside(const TV& location,const ZERO thickness_over_two) const
    {return Lazy_Outside(location);}

    bool Boundary(const TV& location,const T thickness_over_two) const
    {bool strict_inside=location.All_Greater(min_corner+thickness_over_two) && location.All_Less(max_corner-thickness_over_two);
    return !strict_inside && !Outside(location,thickness_over_two);}

    TV Clamp(const TV& location) const
    {return clamp(location,min_corner,max_corner);}

    T Clamp(const T& location) const
    {STATIC_ASSERT(d==1);return Clamp(TV(location)).x;}

    void Enlarge_By_Sign(const TV& v)
    {for(int i=1;i<=d;i++) if(v(i)>0) max_corner(i)+=v(i);else min_corner(i)+=v(i);}

    TV Point_From_Normalized_Coordinates(const TV& weights) const
    {return min_corner+weights*(max_corner-min_corner);}

    bool Contains(const RANGE<TV>& box) const
    {return min_corner.All_Less_Equal(box.min_corner) && max_corner.All_Greater_Equal(box.max_corner);}

    bool Lazy_Intersection(const RANGE<TV>& box) const
    {return min_corner.All_Less_Equal(box.max_corner) && max_corner.All_Greater_Equal(box.min_corner);}

    bool Intersection(const RANGE<TV>& box,const T thickness_over_two) const
    {return Thickened(thickness_over_two).Lazy_Intersection(box);}

    bool Intersection(const RANGE<TV>& box,const ZERO thickness_over_two) const
    {return Lazy_Intersection(box);}

    bool Intersection(const RANGE<TV>& box) const
    {return Lazy_Intersection(box);}

    T Intersection_Area(const RANGE<TV>& box) const
    {return Intersect(*this,box).Robust_Size();}

    void Project_Points_Onto_Line(const TV& direction,T& line_min,T& line_max) const
    {line_min=line_max=TV::Dot_Product(direction,min_corner);TV e=direction*(max_corner-min_corner);
    for(int i=1;i<=d;i++) if(e(i)>0) line_max+=e(i);else line_min+=e(i);}

    static RANGE<TV> Bounding_Box(const TV& p1,const TV& p2)
    {RANGE<TV> box(p1);box.Enlarge_Nonempty_Box_To_Include_Point(p2);return box;}

    static RANGE<TV> Bounding_Box(const TV& p1,const TV& p2,const TV& p3)
    {RANGE<TV> box(p1);box.Enlarge_Nonempty_Box_To_Include_Points(p2,p3);return box;}

    static RANGE<TV> Bounding_Box(const TV& p1,const TV& p2,const TV& p3,const TV& p4)
    {RANGE<TV> box(p1);box.Enlarge_Nonempty_Box_To_Include_Points(p2,p3,p4);return box;}

    template<class T_ARRAY>
    static RANGE<TV> Bounding_Box(const T_ARRAY& points)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    if(!points.Size()) return Empty_Box();
    RANGE<TV> box(points(1));for(int i=2;i<=points.Size();i++) box.Enlarge_Nonempty_Box_To_Include_Point(points(i));return box;}

    RANGE<VECTOR<T,d-1> > Get_Horizontal_Box() const
    {return RANGE<VECTOR<T,d-1> >(min_corner.Horizontal_Vector(),max_corner.Horizontal_Vector());}

    RANGE<VECTOR<T,d-1> > Get_Vertical_Box() const
    {STATIC_ASSERT(d==2);return RANGE<VECTOR<T,d-1> >(min_corner.Vertical_Vector(),max_corner.Vertical_Vector());}

    RANGE<VECTOR<T,d-1> > Remove_Dimension(int dimension) const
    {return RANGE<VECTOR<T,d-1> >(min_corner.Remove_Index(dimension),max_corner.Remove_Index(dimension));}

//#####################################################################
};
template<class TV>
inline RANGE<TV> operator+(const TV& a,const RANGE<TV>& b)
{return RANGE<TV>(a+b.min_corner,a+b.max_corner);}

template<class TV>
inline RANGE<TV> operator-(const TV& a,const RANGE<TV>& b)
{return RANGE<TV>(a-b.max_corner,a-b.min_corner);}

template<class TV>
inline RANGE<TV> operator-(const RANGE<TV>& b,const TV& a)
{return RANGE<TV>(b.min_corner-a,b.max_corner-a);}

template<class TV> inline RANGE<TV> operator*(const typename TV::SCALAR a,const RANGE<TV>& box)
{return box*a;}
}
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#endif
